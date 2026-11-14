#author: haim krupkin+claude
#date: 11.10.2026
#description: step 7b. Stricter OneK1K check of the calibrated cell-type eGenes (step 7 checks only gene overlap).
# For each eGene (calibrated, as step 7: empirical null from the pooled shuffles, BH < 5%) that is also an eQTL
# (q < 0.05) in a matched OneK1K cell type: is it the SAME signal?
#   LD r2 between our lead variant and OneK1K's lead (GTEx 670 donors), and the same direction of effect.
# OneK1K positions are hg19 -> lifted to hg38 (liftOver). OneK1K IDs carry no alleles: orientation by allele frequency
# (OneK1K af vs GTEx ALT frequency at that position; unsure when either is within 0.1 of 0.5 -> direction NA).
# Neutrophils (not in OneK1K, a PBMC study) are compared with classical monocytes.
# Env: CHRS (default all 22), SCAN_TAG (step-6 OUT_TAG, default "").
# Inputs : output/06_scan_<chr>_{genes,top_variants}<SCAN_TAG>.csv, output/genotypes/<chr>, external/onek1k_top_eqtl
# Outputs: output/07b_onek1k_lead_pairs<SCAN_TAG>.csv (per eGene x OneK1K type), 07b_onek1k_lead_summary<SCAN_TAG>.csv
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
onek1k_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/onek1k_top_eqtl/unzipped/OneK1K_TensorQTL_top_eQTL_summary"
LIFTOVER <- "/home/hkrupkin/envs/liftover_env/bin/liftOver"
CHAIN    <- "/labs/smontgom/grps_smontgom/shared/liftOver/chains/hg19ToHg38.over.chain.gz"
CHRS     <- strsplit(Sys.getenv("CHRS", paste0("chr", 1:22, collapse = " ")), " ")[[1]]
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")   # 11.13: main result = 20 hidden-noise PCs; "" = the model without
FDR <- 0.05
ONEK1K_MATCH <- list(Bcells = c("B_IN", "B_MEM", "Plasma"), CD4 = c("CD4_NC", "CD4_ET", "CD4_SOX4"),
                     CD8 = c("CD8_NC", "CD8_ET", "CD8_S100B"), NK = c("NK", "NK_R"),
                     Monocytes = c("Mono_C", "Mono_NC"), Neutrophils = c("Mono_C"))
CELLTYPES <- names(ONEK1K_MATCH)
tmp_dir <- file.path(out_dir, "07b_tmp")
dir.create(tmp_dir, showWarnings = FALSE)

## ================= calibrated eGenes (as step 7) and their lead variant =================
genes <- bind_rows(lapply(CHRS, function(ch) read.csv(file.path(out_dir, paste0("06_scan_", ch, "_genes", SCAN_TAG, ".csv")))))
top   <- bind_rows(lapply(CHRS, function(ch) read.csv(file.path(out_dir, paste0("06_scan_", ch, "_top_variants", SCAN_TAG, ".csv"))))) %>%
  filter(set == "genotype")
hits <- bind_rows(lapply(CELLTYPES, function(ct) {
  d <- genes[genes$test == ct, ]
  nul <- sort(d$p_gene[grepl("^shuffled", d$set)])
  r <- d[d$set == "genotype", ]
  r$fdr_emp <- p.adjust((1 + findInterval(r$p_gene, nul)) / (1 + length(nul)), "BH")
  h <- r$gene_name[r$fdr_emp < FDR]
  t <- top %>% filter(test == ct, gene_name %in% h)
  if (nrow(t) == 0) return(NULL)
  data.frame(chr = t$chr, celltype = ct, gene_name = t$gene_name, our_variant = t$variant_id, our_b = t[[paste0("b_", ct)]])
}))
print(table(hits$celltype))

## ================= OneK1K leads of those genes in the matched types (q < 0.05) =================
okl <- bind_rows(lapply(split(hits, hits$chr), function(hc) {
  ch <- hc$chr[1]
  bind_rows(lapply(unique(unlist(ONEK1K_MATCH[unique(hc$celltype)])), function(ty) {
    f <- file.path(onek1k_dir, paste0("OneK1K_", ty, ".sig_cis_qtl_pairs.", ch, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- fread(f, select = c("phenotype_id", "variant_id", "af", "slope", "qval"))[qval < FDR]
    hc_ty <- hc[sapply(hc$celltype, function(ct) ty %in% ONEK1K_MATCH[[ct]]), ]
    m <- merge(hc_ty, d, by.x = "gene_name", by.y = "phenotype_id")
    if (nrow(m) == 0) return(NULL)
    m$onek1k_type <- ty
    m
  }))
}))
names(okl)[names(okl) == "variant_id"] <- "ok_variant_hg19"
cat("eGene x OneK1K-type pairs:", nrow(okl), "\n")

## ================= liftOver hg19 -> hg38 =================
key <- paste0("k", seq_len(nrow(okl)))
pos19 <- as.integer(sub(".*:", "", okl$ok_variant_hg19))
bed <- data.frame(paste0("chr", sub(":.*", "", okl$ok_variant_hg19)), pos19 - 1L, pos19, key)
fwrite(bed, file.path(tmp_dir, "ok19.bed"), sep = "\t", col.names = FALSE)
system2(LIFTOVER, c(file.path(tmp_dir, "ok19.bed"), CHAIN, file.path(tmp_dir, "ok38.bed"), file.path(tmp_dir, "unmapped.bed")))
b38 <- fread(file.path(tmp_dir, "ok38.bed"), header = FALSE)
okl$pos38 <- b38$V3[match(key, b38$V4)]
cat("lifted to hg38:", sum(!is.na(okl$pos38)), "of", nrow(okl), "\n")

## ================= LD and allele frequency from GTEx genotypes (ALT dosage), one chromosome at a time =================
okl$r <- NA_real_
okl$gtex_alt_af <- NA_real_
for (ch in unique(okl$chr)) {
  pre <- file.path(out_dir, "genotypes", ch)
  pl <- read.plink(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"))
  G <- as(pl$genotypes, "numeric")
  alt <- sapply(strsplit(colnames(G), "_"), `[`, 4)
  flip <- pl$map$allele.2 != alt
  G[, flip] <- 2 - G[, flip]
  for (k in which(okl$chr == ch & !is.na(okl$pos38))) {
    a <- match(okl$our_variant[k], colnames(G))
    b <- which(pl$map$position == okl$pos38[k])
    if (is.na(a) || length(b) != 1) next                    # absent or multi-allelic: skip
    ok_i <- !is.na(G[, a]) & !is.na(G[, b])
    okl$r[k] <- cor(G[ok_i, a], G[ok_i, b])
    okl$gtex_alt_af[k] <- mean(G[, b], na.rm = TRUE) / 2
  }
  rm(G); gc()
}
okl$r2 <- okl$r^2
okl$orient <- with(okl, ifelse(is.na(gtex_alt_af) | abs(af - 0.5) <= 0.1, NA,
                     ifelse(abs(af - gtex_alt_af) < 0.1, 1, ifelse(abs(af - (1 - gtex_alt_af)) < 0.1, -1, NA))))
okl$same_direction <- sign(okl$our_b) == sign(okl$slope * okl$orient * sign(okl$r))

## ================= per eGene: the best-matching OneK1K type; summary per cell type =================
best <- okl %>% group_by(celltype, gene_name) %>% arrange(desc(r2)) %>% slice(1) %>% ungroup()
summary_tab <- best %>% group_by(celltype) %>%
  summarise(egenes_with_onek1k_eqtl = n(),
            testable_ld = sum(!is.na(r2)),
            same_signal_r2_ge_0.8 = sum(r2 >= 0.8, na.rm = TRUE),
            shared_r2_ge_0.5 = sum(r2 >= 0.5, na.rm = TRUE),
            direction_known_r2_ge_0.5 = sum(!is.na(same_direction) & r2 >= 0.5, na.rm = TRUE),
            same_direction_r2_ge_0.5 = sum(same_direction & r2 >= 0.5, na.rm = TRUE),
            .groups = "drop") %>%
  right_join(hits %>% count(celltype, name = "egenes"), by = "celltype") %>%
  relocate(egenes, .after = celltype)
print(as.data.frame(summary_tab))

fwrite(okl, file.path(out_dir, paste0("07b_onek1k_lead_pairs", SCAN_TAG, ".csv")))
fwrite(summary_tab, file.path(out_dir, paste0("07b_onek1k_lead_summary", SCAN_TAG, ".csv")))
unlink(tmp_dir, recursive = TRUE)
message("07b_onek1k_lead_check done")
