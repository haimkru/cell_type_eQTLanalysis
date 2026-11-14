#author: haim krupkin+claude
#date: 11.08.2026
#description: step 7. All results and checks of the cell-type eQTL scan (step 6), over the chromosomes run.
# Tests (per gene; BH < 5% = eGene):
#   <celltype> : the cell type has a genotype effect          het : the 6 cell-type effects differ (cell-type specific)
#   nl         : neutrophil effect differs from the mean lymphoid effect
# eGene calls are CALIBRATED: each gene's p is compared to the pooled gene-level p-values of the 5 shuffled-genotype
# scans (empirical p), then BH. (The model's own p-values were too small for the cell-type-specific test.)
# Checks:
#   A. calibration : each shuffle scored against the other shuffles like real data -> false eGenes (should be ~0)
#   B. replication : the main eGenes, at the same variant, fitted in pre- and post-mortem draws separately:
#                    share with the same sign in both, and the rank correlation of the effects (pre vs post)
#   C. OneK1K      : are the main eGenes of a cell type more often OneK1K eGenes (matching cell type) than the other
#                    tested genes? (gene symbol match; OneK1K has no neutrophils -> classical monocytes as nearest)
#   D. specific genes: for the het eGenes, which cell type has the largest effect (z) at the best variant
#   E. robustness  : the main eGenes refitted with proportions from CIBERSORTx on the RAW mixture (step 6 rawprops):
#                    same sign, rank correlation, share still significant
# Inputs : output/06_scan_<chr>_{genes,top_variants}{,_pre,_post}.csv; external/onek1k_top_eqtl
# Outputs: output/07_egenes.csv, 07_calibration.csv, 07_replication.csv, 07_onek1k.csv, 07_specific_genes.csv
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
onek1k_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/onek1k_top_eqtl/unzipped/OneK1K_TensorQTL_top_eQTL_summary"
CHRS <- strsplit(Sys.getenv("CHRS", "chr21 chr22"), " ")[[1]]
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
TESTS <- c(CELLTYPES, "het", "nl")
FDR <- 0.05
ONEK1K_MATCH <- list(Bcells = c("B_IN", "B_MEM", "Plasma"), CD4 = c("CD4_NC", "CD4_ET", "CD4_SOX4"),
                     CD8 = c("CD8_NC", "CD8_ET", "CD8_S100B"), NK = c("NK", "NK_R"),
                     Monocytes = c("Mono_C", "Mono_NC"), Neutrophils = c("Mono_C"))

## ================= read =================
TEST_TAG <- Sys.getenv("TEST_TAG", "")                       # "_TEST": read the quick-test files of step 6
# 11.10: SCAN_TAG = the OUT_TAG of the step-6 runs to summarise (e.g. "_phi_win100k_shuf25"); also added to the outputs
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")   # 11.13: main result = 20 hidden-noise PCs; "" = the model without
read_scan <- function(kind, tag) {
  bind_rows(lapply(CHRS, function(ch) read.csv(file.path(out_dir, paste0("06_scan_", ch, "_", kind, tag, SCAN_TAG, TEST_TAG, ".csv")))))
}
main <- read_scan("top_variants", "")
pre  <- read_scan("top_variants", "_pre")
post <- read_scan("top_variants", "_post")
rawp <- read_scan("top_variants", "_rawprops")

## ================= gene-level p: Bonferroni (default) or M_eff (11.10) =================
# GENE_P = bonferroni: p_gene = best p x number of cis variants (as written by step 6)
# GENE_P = meff      : p_gene = best p x M_eff, the effective number of independent tests (step 6b, LD-aware);
#                      the same M_eff for the gene's real and shuffled scans, so the calibration below stays valid.
#                      Outputs get the suffix "_meff" (the Bonferroni outputs are kept).
GENE_P <- Sys.getenv("GENE_P", "bonferroni")
stopifnot(GENE_P %in% c("bonferroni", "meff"))
OUT_TAG <- paste0(SCAN_TAG, if (GENE_P == "meff") "_meff" else "")
if (GENE_P == "meff") {
  meff <- bind_rows(lapply(CHRS, function(ch) read.csv(file.path(out_dir, paste0("06b_meff_", ch, ".csv"))))) %>%
    select(gene_name, m_eff)
  regene <- function(d) {
    d <- d %>% select(-any_of("m_eff")) %>% left_join(meff, by = "gene_name")
    stopifnot(!anyNA(d$m_eff))
    d$p_gene <- pmin(d$p_min * d$m_eff, 1)
    d %>% group_by(set, test) %>% mutate(fdr = p.adjust(p_gene, "BH")) %>% ungroup() %>% as.data.frame()
  }
  main <- regene(main)
  cat("gene-level p: best p x M_eff (median M_eff", median(meff$m_eff), ")\n")
}

## ================= empirical p-values against the pooled shuffles =================
# per test: the null = the gene-level p-values (p_gene) of ALL shuffled scans, all chromosomes pooled.
# empirical p of a gene = (1 + number of null values <= its p_gene) / (1 + number of null values); then BH.
# This calibrates the calls even if the model's own p-values are too small for some test.
empirical_fdr <- function(p_obs, p_null) {
  p_null <- sort(p_null)
  n_le <- findInterval(p_obs, p_null)                       # number of null values <= p_obs
  p_emp <- (1 + n_le) / (1 + length(p_null))
  p.adjust(p_emp, "BH")
}
shuffle_sets <- sort(unique(main$set[grepl("^shuffled_", main$set)]))
length(shuffle_sets)
main$fdr_emp <- NA_real_
for (tt in TESTS) {
  null_p <- main$p_gene[main$test == tt & main$set %in% shuffle_sets]
  is_real <- main$test == tt & main$set == "genotype"
  main$fdr_emp[is_real] <- empirical_fdr(main$p_gene[is_real], null_p)
}

## ================= main eGenes =================
egenes <- main %>%
  filter(set == "genotype") %>%
  group_by(test) %>%
  summarise(genes = n(),
            egenes = sum(fdr_emp < FDR),                       # calibrated (empirical null)
            egenes_theoretical = sum(fdr < FDR),               # the model's own p-values, for comparison
            .groups = "drop")
egenes <- egenes[match(TESTS, egenes$test), ]
print(as.data.frame(egenes))

## ================= A. calibration: each shuffle scored against the others =================
# a held-out shuffle treated like real data should give ~0 eGenes
calibration <- bind_rows(lapply(TESTS, function(tt) {
  held_out <- sapply(shuffle_sets, function(s) {
    p_obs <- main$p_gene[main$test == tt & main$set == s]
    p_null <- main$p_gene[main$test == tt & main$set %in% setdiff(shuffle_sets, s)]
    sum(empirical_fdr(p_obs, p_null) < FDR)
  })
  theoretical <- sapply(shuffle_sets, function(s) sum(main$fdr[main$test == tt & main$set == s] < FDR))
  data.frame(test = tt,
             shuffles = length(shuffle_sets),
             false_egenes_calibrated_mean = mean(held_out),
             false_egenes_calibrated_max = max(held_out),
             false_egenes_theoretical_mean = mean(theoretical))
}))
print(calibration %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
main$fdr <- main$fdr_emp                                     # from here on: eGene = calibrated FDR < 5%

## ================= B. replication pre vs post =================
# effect of the tested cell type at the main best variant; for het / nl the neutrophil effect is not meaningful,
# so replication is reported for the 6 cell-type tests only
effect_at <- function(d, ct) {
  d %>% filter(set == "genotype") %>%
    select(gene_name, variant_id, effect = all_of(paste0("b_", ct)), p = all_of(paste0("p_", ct))) %>%
    distinct(gene_name, variant_id, .keep_all = TRUE)
}
replication <- bind_rows(lapply(CELLTYPES, function(ct) {
  hits <- main %>% filter(set == "genotype", test == ct, fdr < FDR) %>% select(gene_name, variant_id)
  if (nrow(hits) == 0) return(data.frame(celltype = ct, main_egenes = 0))
  d <- hits %>%
    left_join(effect_at(pre, ct), by = c("gene_name", "variant_id")) %>%
    rename(effect_pre = effect, p_pre = p) %>%
    left_join(effect_at(post, ct), by = c("gene_name", "variant_id")) %>%
    rename(effect_post = effect, p_post = p) %>%
    filter(!is.na(effect_pre), !is.na(effect_post))
  data.frame(celltype = ct, main_egenes = nrow(hits), with_both_fits = nrow(d),
             same_sign = mean(sign(d$effect_pre) == sign(d$effect_post)),
             rho_pre_vs_post = if (nrow(d) >= 5) cor(d$effect_pre, d$effect_post, method = "spearman") else NA,
             # 11.10: na.rm -- 0.3% of the pre-mortem fits have no robust SE (smaller subset), which made the share NA
             nominal_p05_in_pre = mean(d$p_pre < 0.05, na.rm = TRUE), nominal_p05_in_post = mean(d$p_post < 0.05, na.rm = TRUE),
             n_with_p_pre = sum(!is.na(d$p_pre)), n_with_p_post = sum(!is.na(d$p_post)))
}))
print(replication %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())

## ================= E. robustness: proportions from the raw (not batch-removed) mixture =================
# the main eGenes refitted at the same variant with the other proportions: same sign? still nominally significant?
robustness <- bind_rows(lapply(CELLTYPES, function(ct) {
  hits <- main %>% filter(set == "genotype", test == ct, fdr < FDR) %>%
    select(gene_name, variant_id, effect_main = all_of(paste0("b_", ct)))
  hits <- hits %>% distinct(gene_name, variant_id, .keep_all = TRUE)
  if (nrow(hits) == 0) return(data.frame(celltype = ct, main_egenes = 0))
  d <- hits %>%
    left_join(effect_at(rawp, ct), by = c("gene_name", "variant_id")) %>%
    filter(!is.na(effect))
  data.frame(celltype = ct, main_egenes = nrow(hits), refitted = nrow(d),
             same_sign = mean(sign(d$effect_main) == sign(d$effect)),
             rho_main_vs_raw = if (nrow(d) >= 5) cor(d$effect_main, d$effect, method = "spearman") else NA,
             share_p05_with_raw_props = mean(d$p < 0.05),
             share_p_bonf_with_raw_props = mean(d$p < 0.05 / nrow(d)))
}))
print(robustness %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())

## ================= C. OneK1K =================
read_onek1k <- function(ct, ch) {
  f <- file.path(onek1k_dir, paste0("OneK1K_", ct, ".sig_cis_qtl_pairs.", ch, ".csv"))
  d <- fread(f, select = c("phenotype_id", "qval"))
  d$onek1k_type <- ct
  d
}
onek1k <- bind_rows(lapply(unique(unlist(ONEK1K_MATCH)), function(ct) bind_rows(lapply(CHRS, function(ch) read_onek1k(ct, ch)))))
onek1k_egenes <- onek1k %>% filter(qval < FDR)
onek1k_tab <- bind_rows(lapply(CELLTYPES, function(ct) {
  eg <- unique(onek1k_egenes$phenotype_id[onek1k_egenes$onek1k_type %in% ONEK1K_MATCH[[ct]]])
  d <- main %>% filter(set == "genotype", test == ct)
  d$onek1k <- d$gene_name %in% eg
  d$hit <- d$fdr < FDR
  data.frame(celltype = ct, onek1k_types = paste(ONEK1K_MATCH[[ct]], collapse = "+"),
             our_egenes = sum(d$hit),
             share_egenes_in_onek1k = if (any(d$hit)) mean(d$onek1k[d$hit]) else NA,
             share_other_genes_in_onek1k = mean(d$onek1k[!d$hit]),
             fisher_p = if (any(d$hit)) fisher.test(table(d$hit, d$onek1k))$p.value else NA)
}))
print(onek1k_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())

## ================= D. cell-type-specific genes =================
specific <- main %>% filter(set == "genotype", test == "het", fdr < FDR)
z <- as.matrix(specific[, paste0("b_", CELLTYPES)]) / as.matrix(specific[, paste0("se_", CELLTYPES)])
specific$largest_effect_celltype <- if (nrow(specific) > 0) CELLTYPES[apply(abs(z), 1, which.max)] else character(0)
specific$largest_z <- if (nrow(specific) > 0) apply(abs(z), 1, max, na.rm = TRUE) else numeric(0)
# row maximum ignoring NA (PHI: gated-out cell types are NA); 0 when no other cell type was tested
row_max0 <- function(m) apply(m, 1, function(r) if (all(is.na(r))) 0 else max(r, na.rm = TRUE))
specific <- specific[, c("chr", "gene_name", "variant_id", "p_min", "fdr", "largest_effect_celltype", "largest_z",
                         paste0("b_", CELLTYPES), paste0("p_", CELLTYPES))]
print(table(specific$largest_effect_celltype))

## ================= D2. pattern of each cell type's eGenes (11.10) =================
# the het test has no calibrated genes genome-wide, so per cell type: at each eGene's best variant, is the effect
# concentrated in that cell type? (NOMINAL per-gene p-values; a shared blood-wide eQTL would give ~5% here)
pattern <- bind_rows(lapply(CELLTYPES, function(ct) {
  h <- main %>% filter(set == "genotype", test == ct, fdr < FDR)
  if (nrow(h) == 0) return(data.frame(celltype = ct, egenes = 0))
  others <- setdiff(CELLTYPES, ct)
  z_ct <- h[[paste0("b_", ct)]] / h[[paste0("se_", ct)]]
  z_oth <- as.matrix(h[, paste0("b_", others)]) / as.matrix(h[, paste0("se_", others)])
  data.frame(celltype = ct, egenes = nrow(h),
             median_abs_z_celltype = median(abs(z_ct)),
             median_max_abs_z_others = median(row_max0(abs(z_oth))),
             share_no_other_abs_z_gt2 = mean(row_max0(abs(z_oth)) < 2),
             share_het_p05 = mean(h$p_het < 0.05, na.rm = TRUE),
             # neutrophil-vs-lymphoid contrast: only meaningful for Neutrophils
             share_nl_p05_neut_larger = if (ct == "Neutrophils") mean(h$p_nl < 0.05 & sign(h$nl_est) == sign(h$b_Neutrophils), na.rm = TRUE) else NA)
}))
print(pattern %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())

## ================= save =================
# output file name, with "_meff" when GENE_P = meff
out_csv <- function(name) file.path(out_dir, paste0(name, OUT_TAG, ".csv"))
write.csv(pattern, out_csv("07_celltype_pattern"), row.names = FALSE)
write.csv(egenes, out_csv("07_egenes"), row.names = FALSE)
write.csv(calibration, out_csv("07_calibration"), row.names = FALSE)
write.csv(replication, out_csv("07_replication"), row.names = FALSE)
write.csv(robustness, out_csv("07_robustness_raw_proportions"), row.names = FALSE)
write.csv(onek1k_tab, out_csv("07_onek1k"), row.names = FALSE)
write.csv(specific, out_csv("07_specific_genes"), row.names = FALSE)
out_files <- sapply(c("07_egenes", "07_calibration", "07_replication", "07_onek1k", "07_specific_genes"), out_csv)
out_files
stopifnot(all(file.exists(out_files)))
message("07_results_summary done")
