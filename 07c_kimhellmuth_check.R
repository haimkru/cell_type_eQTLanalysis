#author: haim krupkin+claude
#date: 11.12.2026
#description: step 7c. Our neutrophil eGenes vs the GTEx whole-blood neutrophil interaction eQTLs of Kim-Hellmuth et al.
# 2020 (Science): same 670 donors, different method (one cell type at a time: expression ~ genotype + xCell neutrophil
# score + genotype x score + PEER + genotype PCs; +-1 Mb; eigenMT + BH). Their file: one row per gene (lead variant).
# Checks:
#   1. the two neutrophil measures: our CIBERSORTx neutrophil proportion vs their xCell neutrophil score
#   2. gene overlap: are our neutrophil eGenes more often their neutrophil ieGenes (BH <= 0.05) than our other tested
#      genes? (Fisher; genes tested in both)
#   3. direction: their interaction (b_gi > 0 = the ALT effect grows with neutrophils) vs our neutrophil-minus-lymphoid
#      contrast (nl_est > 0 = the same), at the two lead variants; both use GTEx genotypes with ALT coding, so the
#      orientation is exact (the sign of the LD correlation r between the two leads); r2 from GTEx genotypes
# Inputs : output/06_scan_<chr>_{genes,top_variants}.csv, output/genotypes/<chr>, 05_proportions_6classes.csv,
#          external/kim_hellmuth_ieqtl/ (public GTEx portal files: Whole_Blood.Neutrophils.ieQTL.eigenMT.annotated,
#          xCell scores)
# Outputs: output/07c_kimhellmuth_pairs.csv (per gene), output/07c_kimhellmuth_summary.csv
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
kh_dir   <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/kim_hellmuth_ieqtl"
CHRS <- paste0("chr", 1:22)
FDR  <- 0.05
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")   # step-6 OUT_TAG to check, added to the outputs. 11.13: main result = 20 hidden-noise PCs; "" = the model without

## ================= 1. neutrophil proportion vs xCell neutrophil score =================
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
xcell <- fread(file.path(kh_dir, "GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz"))
setnames(xcell, 1, "cell_type")
xn <- unlist(xcell[cell_type == "Neutrophils", intersect(W_df$SAMPID, names(xcell)), with = FALSE])
both <- intersect(W_df$SAMPID, names(xn))
rho_scores <- cor(W_df$Neutrophils[match(both, W_df$SAMPID)], xn[both], method = "spearman")
cat("samples in both:", length(both), " Spearman rho, our neutrophil proportion vs xCell score:", round(rho_scores, 2), "\n")

## ================= our neutrophil eGenes (calibrated, as step 7) =================
empirical_fdr <- function(p_obs, p_null) {
  p_null <- sort(p_null)
  p.adjust((1 + findInterval(p_obs, p_null)) / (1 + length(p_null)), "BH")
}
genes <- bind_rows(lapply(CHRS, function(ch) fread(file.path(out_dir, paste0("06_scan_", ch, "_genes", SCAN_TAG, ".csv"))) %>%
  filter(test == "Neutrophils")))
ours <- genes %>% filter(set == "genotype")
ours$fdr_emp <- empirical_fdr(ours$p_gene, genes$p_gene[genes$set != "genotype"])
top <- bind_rows(lapply(CHRS, function(ch) fread(file.path(out_dir, paste0("06_scan_", ch, "_top_variants", SCAN_TAG, ".csv"))) %>%
  filter(set == "genotype", test == "Neutrophils") %>% select(chr, gene_name, variant_id, b_Neutrophils, nl_est, p_nl)))
ours <- ours %>% select(gene_name, fdr_emp) %>% left_join(top, by = "gene_name")

## ================= their neutrophil ieGenes =================
kh <- fread(file.path(kh_dir, "GTEx_Analysis_v8_ieQTL", "Whole_Blood.Neutrophils.ieQTL.eigenMT.annotated.txt.gz")) %>%
  select(gene_name, kh_variant = variant_id, biotype, b_gi, pval_gi, pval_adj_bh) %>%
  group_by(gene_name) %>% filter(n() == 1) %>% ungroup()           # unique symbols only
cat("their ieGenes (BH <= 0.05), all biotypes:", sum(kh$pval_adj_bh <= FDR),
    "; protein coding + lincRNA:", sum(kh$pval_adj_bh <= FDR & kh$biotype %in% c("protein_coding", "lincRNA")), "\n")

## ================= 2. gene overlap =================
m <- inner_join(ours, kh, by = "gene_name") %>%
  mutate(our_egene = fdr_emp < FDR, their_iegene = pval_adj_bh <= FDR)
tab <- table(ours = m$our_egene, theirs = m$their_iegene)
print(tab)
ft <- fisher.test(tab)
overlap <- data.frame(genes_in_both = nrow(m),
                      our_egenes = sum(m$our_egene), their_iegenes = sum(m$their_iegene),
                      both = sum(m$our_egene & m$their_iegene),
                      share_of_ours_in_theirs = mean(m$their_iegene[m$our_egene]),
                      share_of_other_genes_in_theirs = mean(m$their_iegene[!m$our_egene]),
                      share_of_theirs_in_ours = mean(m$our_egene[m$their_iegene]),
                      odds_ratio = unname(ft$estimate), fisher_p = ft$p.value)
print(overlap, digits = 3)

## ================= 3. direction at the two lead variants (genes called by both) =================
pairs <- m %>% filter(our_egene, their_iegene)
pairs$r <- NA_real_
for (ch in unique(pairs$chr)) {
  pre <- file.path(out_dir, "genotypes", ch)
  bim <- fread(paste0(pre, ".bim"), header = FALSE)
  rows <- which(pairs$chr == ch)
  need <- unique(c(pairs$variant_id[rows], pairs$kh_variant[rows]))
  k <- which(bim$V2 %in% need)
  if (length(k) == 0) next
  pl <- read.plink(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"), select.snps = k)
  G <- as(pl$genotypes, "numeric")
  alt <- sapply(strsplit(colnames(G), "_"), `[`, 4)
  flip <- pl$map$allele.2 != alt
  G[, flip] <- 2 - G[, flip]                                          # ALT dosage, as step 6
  for (i in rows) {
    if (pairs$variant_id[i] %in% colnames(G) && pairs$kh_variant[i] %in% colnames(G)) {
      pairs$r[i] <- cor(G[, pairs$variant_id[i]], G[, pairs$kh_variant[i]], use = "complete.obs")
    }
  }
}
pairs <- pairs %>%
  mutate(r2 = r^2,
         same_direction = sign(nl_est) == sign(b_gi * sign(r)))      # their effect oriented to our lead variant
direction <- pairs %>%
  summarise(genes_called_by_both = n(),
            testable_ld = sum(!is.na(r)),
            same_lead_r2_ge_0.8 = sum(r2 >= 0.8, na.rm = TRUE),
            shared_r2_ge_0.5 = sum(r2 >= 0.5, na.rm = TRUE),
            same_direction_r2_ge_0.5 = sum(same_direction & r2 >= 0.5, na.rm = TRUE),
            same_direction_r2_ge_0.5_and_our_nl_p05 = sum(same_direction & r2 >= 0.5 & p_nl < 0.05, na.rm = TRUE),
            shared_and_our_nl_p05 = sum(r2 >= 0.5 & p_nl < 0.05, na.rm = TRUE))
print(as.data.frame(direction))

summary_tab <- bind_cols(data.frame(rho_neutrophil_proportion_vs_xcell = rho_scores), overlap, direction)
fwrite(pairs, file.path(out_dir, paste0("07c_kimhellmuth_pairs", SCAN_TAG, ".csv")))
fwrite(summary_tab, file.path(out_dir, paste0("07c_kimhellmuth_summary", SCAN_TAG, ".csv")))
cat("outputs:", file.path(out_dir, paste0("07c_kimhellmuth_pairs", SCAN_TAG, ".csv")), file.path(out_dir, paste0("07c_kimhellmuth_summary", SCAN_TAG, ".csv")), "\n")
message("07c_kimhellmuth_check done")
