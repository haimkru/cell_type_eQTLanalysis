#author: haim krupkin+claude
#date: 11.13.2026
#description: step 7e. Variants for the leakage check WITHOUT winner's curse.
# Problem: when a gene and its variant are chosen because they are significant in OUR cell type A, A's effect is
# inflated (winner's curse), and the A-vs-B effect correlation is biased. Fix: choose genes and variants in a THIRD
# dataset. Here: OneK1K (single-cell eQTLs, 982 donors, no overlap with GTEx).
#   gene    = a OneK1K eGene (q < 0.05) in any OneK1K type matched to our 6 cell types
#   variant = that gene's OneK1K lead variant (from the OneK1K type where it is strongest), lifted hg19 -> hg38,
#             matched to our variant by position
# Step 6 is then run with MODEL=lead at these variants (one variant per gene, all 6 effects estimated), and
# 09_results_figures.R (section 16) correlates the effects of two cell types with the r_b correction (Qi et al. 2018).
# Both effects of a pair come from OUR data at the same variant, so allele orientation does not matter.
# Inputs : OneK1K top-eQTL summary files (external/onek1k_top_eqtl), output/genotypes/<chr>.bim
# Outputs: output/07e_onek1k_lead_variants.csv  (gene_name, variant_id, onek1k_type: one row per gene; the LEAD_FILE)
#          output/07e_onek1k_egenes_by_celltype.csv  (gene_name, celltype: in which of our cell types' matched
#          OneK1K types the gene is an eGene; used to pick the genes of each pair)
# Run: Rscript 07e_onek1k_lead_variants.R ; then per chromosome:
#   CHR=chr22 MODEL=lead WINDOW=1e6 LEAD_FILE=output/07e_onek1k_lead_variants.csv OUT_TAG=_epc20_ok1k Rscript 06_cis_scan.R
#   (WINDOW=1e6: OneK1K used a 1 Mb window, so its leads can be up to 1 Mb from the TSS)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
onek1k_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/onek1k_top_eqtl/unzipped/OneK1K_TensorQTL_top_eQTL_summary"
LIFTOVER <- "/home/hkrupkin/envs/liftover_env/bin/liftOver"
CHAIN    <- "/labs/smontgom/grps_smontgom/shared/liftOver/chains/hg19ToHg38.over.chain.gz"
tmp_dir  <- file.path(out_dir, "tmp_07e")
dir.create(tmp_dir, showWarnings = FALSE)
CHRS <- paste0("chr", 1:22)
FDR  <- 0.05
# our cell type -> the OneK1K types it is matched to (as in step 7 / 7b; neutrophils: classical monocytes)
ONEK1K_MATCH <- list(Bcells = c("B_IN", "B_MEM", "Plasma"), CD4 = c("CD4_NC", "CD4_ET", "CD4_SOX4"),
                     CD8 = c("CD8_NC", "CD8_ET", "CD8_S100B"), NK = c("NK", "NK_R"),
                     Monocytes = c("Mono_C", "Mono_NC"), Neutrophils = c("Mono_C"))
match_tab <- data.frame(celltype = rep(names(ONEK1K_MATCH), lengths(ONEK1K_MATCH)),
                        onek1k_type = unlist(ONEK1K_MATCH), row.names = NULL)
onek1k_types <- unique(match_tab$onek1k_type)

## ================= 1. OneK1K eGenes (q < 0.05) with their lead variant, every type and chromosome =================
ok_list <- list()
for (ty in onek1k_types) {
  for (ch in CHRS) {
    f <- file.path(onek1k_dir, paste0("OneK1K_", ty, ".sig_cis_qtl_pairs.", ch, ".csv"))
    if (!file.exists(f)) next
    d <- fread(f, select = c("phenotype_id", "variant_id", "pval_nominal", "qval"))
    d <- d %>% filter(qval < FDR)
    d <- d %>% mutate(onek1k_type = ty, chr = ch)
    ok_list[[paste(ty, ch)]] <- d
  }
}
ok <- bind_rows(ok_list)
ok <- ok %>% rename(gene_name = phenotype_id, variant_hg19 = variant_id)
cat("OneK1K eGene x type rows:", nrow(ok), "; genes:", n_distinct(ok$gene_name), "\n")

## ================= 2. in which of OUR cell types is each gene a (matched) OneK1K eGene? =================
egenes_by_celltype <- ok %>%
  inner_join(match_tab, by = "onek1k_type", relationship = "many-to-many") %>%
  distinct(gene_name, celltype)
print(table(egenes_by_celltype$celltype))
write.csv(egenes_by_celltype, file.path(out_dir, "07e_onek1k_egenes_by_celltype.csv"), row.names = FALSE)

## ================= 3. one variant per gene: the lead of the OneK1K type where the gene is strongest =================
lead <- ok %>%
  group_by(gene_name) %>%
  arrange(pval_nominal) %>%
  slice(1) %>%
  ungroup()
cat("genes with one OneK1K lead variant:", nrow(lead), "\n")

## ================= 4. liftOver hg19 -> hg38 =================
lead <- lead %>% mutate(key = paste0("k", row_number()),
                        pos19 = as.integer(sub(".*:", "", variant_hg19)))
bed <- data.frame(chrom = lead$chr, start = lead$pos19 - 1L, end = lead$pos19, key = lead$key)
fwrite(bed, file.path(tmp_dir, "ok19.bed"), sep = "\t", col.names = FALSE)
system2(LIFTOVER, c(file.path(tmp_dir, "ok19.bed"), CHAIN, file.path(tmp_dir, "ok38.bed"), file.path(tmp_dir, "unmapped.bed")))
b38 <- fread(file.path(tmp_dir, "ok38.bed"), header = FALSE)
lead$chr38 <- b38$V1[match(lead$key, b38$V4)]
lead$pos38 <- b38$V3[match(lead$key, b38$V4)]
cat("lifted to hg38:", sum(!is.na(lead$pos38)), "of", nrow(lead), "\n")
# a lead lifted to another chromosome is dropped
lead <- lead %>% filter(!is.na(pos38), chr38 == chr)

## ================= 5. our variant at that position (GTEx genotypes; a position with 2+ variants is dropped) =================
match_list <- list()
for (ch in unique(lead$chr)) {
  bim <- fread(file.path(out_dir, "genotypes", paste0(ch, ".bim")), header = FALSE, select = c(2, 4))
  names(bim) <- c("variant_id", "pos38")
  bim <- bim %>% group_by(pos38) %>% filter(n() == 1) %>% ungroup()     # one variant per position only
  m <- lead %>% filter(chr == ch) %>% inner_join(bim, by = "pos38")
  match_list[[ch]] <- m
}
lead_matched <- bind_rows(match_list)
cat("matched to a GTEx variant:", nrow(lead_matched), "genes\n")

## ================= save: the LEAD_FILE for step 6 (MODEL=lead) =================
lead_out <- lead_matched %>% select(gene_name, variant_id, onek1k_type, variant_hg19)
write.csv(lead_out, file.path(out_dir, "07e_onek1k_lead_variants.csv"), row.names = FALSE)
cat("saved:", file.path(out_dir, "07e_onek1k_lead_variants.csv"), "\n")
