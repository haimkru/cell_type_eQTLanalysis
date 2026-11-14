#author: haim krupkin+claude
#date: 11.14.2026
#description: step 12b. Calibrated cell-type trans-eQTL calls (step 12) and checks against the usual trans artefacts.
# Calls, per test (6 cell types, het, nl), all chromosomes pooled:
#   empirical FDR at a p threshold t = (mean number of shuffled pairs with p <= t, over the N shuffles)
#                                      / (number of real pairs with p <= t)
#   (scaled by real / mean shuffled number of tests, ~1), made monotone in p. Trans-eQTL = FDR < 5%.
#   The shuffles keep the LD between candidates and the correlation between genes, so this FDR accounts for both.
# Checks:
#   A. calibration : each shuffle scored against the other shuffles like real data -> false calls (should be ~0)
#   B. proportion QTL: a variant that changes a cell-type proportion can create false trans hits in many genes (whatever
#      the proportion estimate misses ends up in the residual). Per called variant: lm(proportion ~ genotype +
#      covariates) for each cell type; flagged if p < 0.05 / number of called variants.
#   C. hubs        : variants with >= HUB_N trans genes in one test (real master regulators look like this, but so do
#                    proportion QTLs) -> listed, with their proportion-QTL p-values.
#   D. cell-type specific: for each cell-type call, are the other cell types' effects at the pair small (|z| < 2)?
#   E. cis link    : for candidates that are cis leads (step 12a), is the trans effect in the same cell type as the cis one?
#   F. cross-mappability (Saha & Battle 2018, default on): pairs whose trans gene is cross-mappable with any gene within
#      1 Mb of the variant are flagged (misaligned reads look like trans).
# Inputs : output/12_trans_<chr>_{pairs,ntests}<SCAN_TAG>.csv, 12a_trans_candidates<SCAN_TAG>.csv, step-6 inputs
# Outputs: output/12b_trans_calls<SCAN_TAG>.csv, 12b_trans_egenes<SCAN_TAG>.csv, 12b_trans_calibration<SCAN_TAG>.csv,
#          12b_trans_variants<SCAN_TAG>.csv (per called variant: trans genes per test, proportion-QTL p, hub flag)
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHRS     <- strsplit(Sys.getenv("CHRS", "chr21 chr22"), " ")[[1]]   # chromosomes of the GENES scanned in step 12
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")
CAND_FILE <- Sys.getenv("CAND_FILE", file.path(out_dir, paste0("12a_trans_candidates", SCAN_TAG, ".csv")))
CROSSMAP_FILE <- Sys.getenv("CROSSMAP_FILE", "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/crossmap_saha_battle_2018/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz")
pc_file  <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
TESTS    <- c(CELLTYPES, "het", "nl")
FDR      <- 0.05
HUB_N    <- 10
OTHER_COV <- c("post_mortem", "ventilator", "SMTSISCH", "SMRIN", "centre_B1", "male", "age")

## ================= read step 12 =================
pairs <- bind_rows(lapply(CHRS, function(ch) fread(file.path(out_dir, paste0("12_trans_", ch, "_pairs", SCAN_TAG, ".csv")))))
ntests <- bind_rows(lapply(CHRS, function(ch) fread(file.path(out_dir, paste0("12_trans_", ch, "_ntests", SCAN_TAG, ".csv")))))
ntests <- ntests %>% group_by(set, test) %>% summarise(n_tests = sum(n_tests), .groups = "drop")
shuffle_sets <- sort(unique(ntests$set[grepl("^shuffled_", ntests$set)]))
cat("genes' chromosomes:", length(CHRS), " shuffles:", length(shuffle_sets), "\n")
print(as.data.frame(ntests %>% filter(set == "genotype")))

## ================= empirical FDR against the shuffles =================
# p_obs: p-values of the "real" set (only the kept ones, all < P_KEEP); p_null: those of the null sets pooled;
# n_null_sets: how many null sets were pooled; ratio: real tests / mean null tests
empirical_fdr_trans <- function(p_obs, p_null, n_null_sets, ratio) {
  o <- order(p_obs)
  p_sorted <- p_obs[o]
  n_real_le <- seq_along(p_sorted)                            # real pairs with p <= each threshold
  n_null_le <- findInterval(p_sorted, sort(p_null))          # null pairs (all null sets) with p <= threshold
  fdr <- (n_null_le / n_null_sets) * ratio / n_real_le
  fdr <- pmin(rev(cummin(rev(fdr))), 1)                       # monotone: never larger than at a bigger threshold
  out <- numeric(length(p_obs))
  out[o] <- fdr
  out
}
p_of <- function(set_names, tt) {
  col <- paste0("p_", tt)
  v <- pairs[[col]][pairs$set %in% set_names]
  v[!is.na(v)]
}
n_of <- function(set_names, tt) mean(ntests$n_tests[ntests$set %in% set_names & ntests$test == tt])

calls <- list()
for (tt in TESTS) {
  col <- paste0("p_", tt)
  real <- pairs %>% filter(set == "genotype", !is.na(.data[[col]]))
  if (nrow(real) == 0) next
  ratio <- n_of("genotype", tt) / n_of(shuffle_sets, tt)
  real$fdr_emp <- empirical_fdr_trans(real[[col]], p_of(shuffle_sets, tt), length(shuffle_sets), ratio)
  real$test <- tt
  real$p_test <- real[[col]]
  calls[[tt]] <- real %>% filter(fdr_emp < FDR)
}
calls <- as.data.frame(bind_rows(calls))
egenes <- bind_rows(lapply(TESTS, function(tt) {
  d <- calls %>% filter(test == tt)
  data.frame(test = tt, tests = n_of("genotype", tt), trans_pairs = nrow(d),
             trans_genes = n_distinct(d$gene_name), trans_variants = n_distinct(d$variant_id),
             largest_p_called = if (nrow(d) > 0) max(d$p_test) else NA)
}))
print(egenes)

## ================= A. calibration: each shuffle as if real =================
calibration <- bind_rows(lapply(TESTS, function(tt) {
  false_calls <- sapply(shuffle_sets, function(s) {
    others <- setdiff(shuffle_sets, s)
    p_obs <- p_of(s, tt)
    if (length(p_obs) == 0) return(0)
    sum(empirical_fdr_trans(p_obs, p_of(others, tt), length(others), n_of(s, tt) / n_of(others, tt)) < FDR)
  })
  data.frame(test = tt, shuffles = length(shuffle_sets),
             false_pairs_mean = mean(false_calls), false_pairs_max = max(false_calls))
}))
print(calibration)

## ================= B. proportion QTL of each called variant =================
called_vars <- unique(calls$variant_id)
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
W <- as.matrix(W_df[match(cov_df$SAMPID, W_df$SAMPID), CELLTYPES])
pcs <- fread(pc_file)
setnames(pcs, 1, "ID")
pcs <- pcs[ID %in% paste0("PC", 1:5)]
PC <- t(as.matrix(pcs[, cov_df$SUBJID, with = FALSE]))
covs <- cbind(as.matrix(cov_df[, OTHER_COV]), PC)
read_geno <- function(ids) {                                 # ALT dosage, samples x variants
  ch <- sub("_.*", "", ids[1])
  pfx <- file.path(out_dir, "genotypes", ch)
  pl <- read.plink(paste0(pfx, ".bed"), paste0(pfx, ".bim"), paste0(pfx, ".fam"), select.snps = ids)
  G <- as(pl$genotypes, "numeric")
  flip <- pl$map$allele.2 != sapply(strsplit(colnames(G), "_"), `[`, 4)
  G[, flip] <- 2 - G[, flip]
  G[match(cov_df$SUBJID, rownames(G)), , drop = FALSE]
}
var_tab <- data.frame(variant_id = called_vars)
if (length(called_vars) > 0) {
  G <- do.call(cbind, lapply(split(called_vars, sub("_.*", "", called_vars)), read_geno))
  # p of the genotype in lm(proportion ~ genotype + covariates), per variant and cell type
  pqtl <- t(sapply(called_vars, function(v) {
    sapply(CELLTYPES, function(ct) {
      fit <- lm(W[, ct] ~ G[, v] + covs)
      summary(fit)$coefficients[2, 4]
    })
  }))
  colnames(pqtl) <- paste0("propqtl_p_", CELLTYPES)
  var_tab <- cbind(var_tab, pqtl)
  var_tab$propqtl_min_p <- apply(pqtl, 1, min)
  var_tab$propqtl_flag <- var_tab$propqtl_min_p < 0.05 / length(called_vars)
  # C. trans genes per variant and test; hub = >= HUB_N in any test
  per_test <- calls %>% count(variant_id, test) %>% pivot_wider(names_from = test, values_from = n, values_fill = 0,
                                                                 names_prefix = "n_trans_")
  var_tab <- var_tab %>% left_join(per_test, by = "variant_id")
  n_cols <- grep("^n_trans_", names(var_tab), value = TRUE)
  var_tab$hub <- apply(var_tab[, n_cols, drop = FALSE], 1, max) >= HUB_N
  cat("called variants:", nrow(var_tab), " proportion-QTL flagged:", sum(var_tab$propqtl_flag),
      " hubs (>=", HUB_N, "genes):", sum(var_tab$hub), "\n")
}

## ================= D-F: only when there are calls =================
if (nrow(calls) == 0) {
  cat("no trans-eQTL calls at FDR <", FDR, "\n")
} else {

## ================= D. cell-type specific: other cell types small at the pair =================
row_max0 <- function(m) apply(m, 1, function(r) if (all(is.na(r))) 0 else max(r, na.rm = TRUE))
calls$max_abs_z_other <- NA_real_
for (ct in CELLTYPES) {
  i <- which(calls$test == ct)
  if (length(i) == 0) next
  others <- setdiff(CELLTYPES, ct)
  z <- as.matrix(calls[i, paste0("b_", others), drop = FALSE]) / as.matrix(calls[i, paste0("se_", others), drop = FALSE])
  calls$max_abs_z_other[i] <- row_max0(abs(z))
}
calls$specific <- calls$test %in% CELLTYPES & calls$max_abs_z_other < 2

## ================= E. cis link (cis-lead candidates) =================
cand <- read.csv(CAND_FILE)
calls <- calls %>% left_join(cand %>% select(variant_id, cis_source = source, cis_genes), by = "variant_id")
calls$same_celltype_as_cis <- mapply(function(tt, src) tt %in% strsplit(src, ";")[[1]], calls$test, calls$cis_source)

## ================= F. cross-mappability (Saha & Battle 2018) =================
# Their recommendation: drop a trans pair if the trans gene is cross-mappable (> 0 shared k-mers, up to 2 mismatches)
# with ANY gene within 1 Mb of the variant: reads from a local gene that mis-map to the distant gene look like trans.
# File: GRCh38 / GENCODE v26 (the annotation GTEx v8 used), symmetric cross-mappability, Ensembl ids; downloaded 11.14
# from figshare 7315136 (md5 checked). CROSSMAP_FILE="" switches the check off.
gtf_file <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf"
gtf <- fread(cmd = paste("grep -P '\\tgene\\t'", gtf_file), header = FALSE)
gtf$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)
gtf$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9)
gtf$gene_type <- sub('.*gene_type "([^"]+)".*', "\\1", gtf$V9)
CROSSMAP_WIN <- 1e6
calls$crossmappable <- NA
if (nzchar(CROSSMAP_FILE)) {
  # Ensembl ids of each called pair's trans gene and of the genes within 1 Mb of its variant
  calls$var_chr <- sub("_.*", "", calls$variant_id)
  calls$var_pos <- as.numeric(sapply(strsplit(calls$variant_id, "_"), `[`, 2))
  trans_ids <- lapply(calls$gene_name, function(g) gtf$gene_id[gtf$gene_name == g])
  local_ids <- lapply(seq_len(nrow(calls)), function(i) {
    gtf$gene_id[gtf$V1 == calls$var_chr[i] & gtf$V5 >= calls$var_pos[i] - CROSSMAP_WIN & gtf$V4 <= calls$var_pos[i] + CROSSMAP_WIN]
  })
  # the cross-mappability file, only the rows touching a called trans gene (28M rows in all)
  cm <- fread(CROSSMAP_FILE, header = FALSE, select = 1:2, col.names = c("gene1", "gene2"))
  needed <- unique(unlist(trans_ids))
  cm <- cm[gene1 %in% needed | gene2 %in% needed]
  cm_key <- unique(c(paste(cm$gene1, cm$gene2), paste(cm$gene2, cm$gene1)))
  rm(cm)
  calls$crossmappable <- mapply(function(tr, loc) {
    any(outer(tr, loc, paste) %in% cm_key)
  }, trans_ids, local_ids)
  cat("calls cross-mappable with a gene within 1 Mb of the variant:", sum(calls$crossmappable), "of", nrow(calls), "\n")
}
calls <- calls %>% left_join(var_tab %>% select(variant_id, propqtl_min_p, propqtl_flag, hub), by = "variant_id")

## ================= G. built-in mapping-artefact flags (11.14) =================
# chr22 test: the strongest call (p 5e-44) was SIRPAP1, a processed pseudogene of SIRPA, at SIRPA's own cis lead:
# reads of SIRPA mis-mapped to its pseudogene. So, without needing a cross-mappability file:
#   pseudogene : the trans gene is a pseudogene (GENCODE gene_type contains "pseudogene")
#   name_link  : the trans gene's name starts with the name of a cis gene of the variant, or the reverse
#                (SIRPA -> SIRPAP1; paralog / pseudogene naming)
#   receptor   : immunoglobulin, T-cell receptor, KIR or HLA genes (hypervariable, hard to map, clonal)
gene_type <- setNames(gtf$gene_type, gtf$gene_name)
calls$gene_type <- unname(gene_type[calls$gene_name])  # (gtf read in F)
calls$pseudogene <- grepl("pseudogene", calls$gene_type)
calls$name_link <- mapply(function(g, cis) {
  if (is.na(cis)) return(FALSE)
  cis <- strsplit(cis, ";")[[1]]
  any(startsWith(g, cis) | startsWith(cis, g))
}, calls$gene_name, calls$cis_genes)
calls$receptor <- grepl("^IG[HKL][VDJC]|^TR[ABDG][VDJC]|^KIR|^HLA-", calls$gene_name)
# clean call = none of the artefact flags (cross-mappability only if the file was given)
calls$clean <- !calls$pseudogene & !calls$name_link & !calls$receptor & !calls$propqtl_flag &
  !(calls$crossmappable %in% TRUE)
cat("calls:", nrow(calls), " clean (no pseudogene / name link / receptor / proportion-QTL / cross-mappable flag):",
    sum(calls$clean), "\n")
}

## ================= summary per test =================
if (nrow(calls) > 0) {
  check_tab <- calls %>%
    group_by(test) %>%
    summarise(trans_pairs = n(),
              share_specific = mean(specific),
              share_same_celltype_as_cis = mean(same_celltype_as_cis, na.rm = TRUE),
              share_variant_propqtl = mean(propqtl_flag),
              share_variant_hub = mean(hub),
              share_crossmappable = mean(crossmappable),
              share_pseudogene = mean(pseudogene),
              share_name_link = mean(name_link),
              share_receptor = mean(receptor),
              clean_pairs = sum(clean),
              clean_genes = n_distinct(gene_name[clean]),
              .groups = "drop")
  print(as.data.frame(check_tab))
  egenes <- egenes %>% left_join(check_tab %>% select(-trans_pairs), by = "test")
}

## ================= save =================
out_csv <- function(name) file.path(out_dir, paste0(name, SCAN_TAG, ".csv"))
write.csv(calls, out_csv("12b_trans_calls"), row.names = FALSE)
write.csv(egenes, out_csv("12b_trans_egenes"), row.names = FALSE)
write.csv(calibration, out_csv("12b_trans_calibration"), row.names = FALSE)
write.csv(var_tab, out_csv("12b_trans_variants"), row.names = FALSE)
sapply(c("12b_trans_calls", "12b_trans_egenes", "12b_trans_calibration", "12b_trans_variants"), out_csv)
message("12b_trans_summary done")
