#author: haim krupkin+claude
#date: 11.07.2026 (11.08: robust HC3 errors; MODEL pre/post; renumbered from step 10 to step 6)
#manually verified by hk this part is hella complex.... 11/14/2026
#this is gonna take some time to verifiy and the code long as hell.


#description: step 6. Cell-type eQTL scan over ALL cis variants of ALL expressed genes of one chromosome (CHR), with
# contrast tests between cell types. ONE joint model per gene and variant: all 6 cell-type effects estimated together.
# Model (TCA, Rahmani 2019), per sample i and cell type h:
#   cell-type expression  Z_h(i) = mu_h + genotype_i * gamma_h + death_i * d_h + noise_h   (noise SD sigma_h)
#   bulk expression       X(i)   = sum_h W_h(i) * Z_h(i) + covariates_i + genotype_i x death_i + noise (SD tau)
#   gamma_h = the genotype effect in cell type h, learned from how the proportions W vary between people.
#   Death (blood draw, ventilator) has its own effect per cell type, and genotype x death is a bulk covariate, so a
#   cell-type effect can only come from proportion differences WITHIN death type.
# MODEL: deathadj (main: all samples, all variants) / pre / post (replication: only pre- or post-mortem draws, only
# the main scan's best variants; death variables that do not vary in the group are dropped) / rawprops (robustness:
# all samples, the main scan's best variants, proportions from CIBERSORTx on the raw, not batch-removed, mixture).
# How it is computed (fast):
#   - TCA estimates the effects by weighted least squares, each sample weighted by its expected variance #thats true
#     sum_h W_h^2 sigma_h^2 + tau^2 (TCA's tca.fit_means_vars). So:
#     a. one TCA fit per gene WITHOUT genotype (death per cell type in C1, other covariates in C2) -> sigma_h, tau
#     b. for every variant: the same weighted regression with the genotype added
#        design = [W_h | W_h x post_mortem | W_h x ventilator | covariates | W_h x genotype | genotype x death]
#     (checked: effects r = 0.995 and -log10 p r = 0.98 vs TCA's own per-gene fit)
#   - standard errors: robust HC3 sandwich (plain errors were too small for CD4 and the contrasts: a few samples with
#     extreme proportions have huge leverage; shuffled genotypes gave 12-18 false eGenes, HC3: 0-1)
#   - the genotype effect of each cell type = the W_h x genotype coefficients; their covariance gives contrasts:
#        heterogeneity : are the 6 cell-type effects all equal? (Wald, 5 df)  <- "cell-type specific"
#        neut_vs_lymph : neutrophil effect minus the mean of B, CD4, CD8, NK (Wald, 1 df)
#   - null: the same scan with the genotype rows shuffled across samples (5 shuffles, LD kept). Step 7 pools them
#     into an empirical null per test, so the final eGene calls are calibrated whatever the model's p-values do.
# Gene level (per test): smallest p over the gene's variants x number of variants (Bonferroni within gene, max 1),
# then BH across genes; eGene = BH < 5%.
# Variants: within 100 kb of the TSS (11.10; was 1 Mb), MAF >= 0.05 (interactions with rarer variants are unstable).
# Genotype = ALT-allele dosage. Expression = log2(TPM + 1), Hb-free (step 1). Proportions = step 5.
# Inputs : output/01_tpm_noHb.rds, 01_gene_subset.txt, 02_covariates.csv, 05_proportions_6classes.csv;
#          output/genotypes/<CHR>.{bed,bim,fam} (made by the master); GENCODE v26; GTEx genotype PCs
# Outputs (<m> = "" for deathadj, "_pre", "_post"):
#          output/06_scan_<CHR>_genes<m>.csv        (gene x test: best variant, n variants, gene-level p, BH; real + shuffled)
#          output/06_scan_<CHR>_top_variants<m>.csv (the same rows with all 6 effects, SEs and p-values)
#          output/06_scan_<CHR><m>_null_fit.rds     (TCA fit without genotype: sigma per gene and cell type, tau)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))
library(TCA)
library(parallel)

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHR      <- Sys.getenv("CHR", "chr22")
# genotypes of this chromosome: extracted by the master (step 6) from the genome-wide GTEx plink file
geno_prefix <- file.path(out_dir, "genotypes", CHR)
gtex     <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8"
gtf_file <- file.path(gtex, "references/gencode.v26.GRCh38.genes.gtf")
pc_file  <- file.path(gtex, "eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt")
attr_ann_file <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/edQTLs_GTEx/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
# cis window around the TSS. 11.10: default 100 kb (was 1 Mb). Tested on chr21-22 with 25 shuffles: neutrophil eGenes
# 21 -> 30, first calibrated het / neutrophil-vs-lymphoid genes (0 -> 2), calibration unchanged (~0.1 false per shuffle).
# The 1 Mb genome-wide results are kept in output/superseded/window1Mb_1110/.
WINDOW   <- as.numeric(Sys.getenv("WINDOW", "1e5"))
MIN_MAF  <- 0.05 #hk true 11/14/2026
N_CORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
N_TEST   <- as.integer(Sys.getenv("N_TEST", "0"))           # > 0: only the first N_TEST genes (quick test)
# deathadj (main, all samples) / pre / post (replication) / rawprops (robustness: proportions from the raw mixture)
MODEL    <- Sys.getenv("MODEL", "deathadj")
# 11.10: spike (power check, step 8): the deathadj model on real expression with KNOWN cell-type effects added
# 11.10: lead (alternative design, under test): the deathadj model at ONE variant per gene, the GTEx v8 bulk
#   whole-blood lead eQTL variant (no x ~4,450 within-gene Bonferroni), still with 5 shuffles for calibration
stopifnot(MODEL %in% c("deathadj", "pre", "post", "rawprops", "spike", "lead"))
model_tag <- if (MODEL == "deathadj") "" else paste0("_", MODEL)
# 11.13: MAIN MODEL = with N_EPC = 20 hidden-noise PCs (step 5c) as covariates (user's decision; chr21-22 and genome-wide
# tests 11.12-13: about twice the eGenes, all checks hold). Its outputs carry OUT_TAG "_epc20" (default when N_EPC = 20).
# The model without them: N_EPC=0 OUT_TAG="" (those 11.11 results are kept, untagged).
N_EPC   <- as.integer(Sys.getenv("N_EPC", "20"))
OUT_TAG <- Sys.getenv("OUT_TAG", if (N_EPC > 0) paste0("_epc", N_EPC) else "")
SPIKE_FILE <- file.path(out_dir, paste0("08_spike_plan_", CHR, OUT_TAG, ".csv"))   # gene_name, celltype, beta_sd (step 8)
# MODEL = lead: one variant per gene from this file (columns gene_name, variant_id). Default: GTEx bulk lead variants.
# 11.13: step 7e writes OneK1K lead variants (selected outside our data, for the leakage check without winner's curse)
LEAD_FILE  <- Sys.getenv("LEAD_FILE", "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz")
# shuffled-genotype null scans: 5 for the main model (pooled into an empirical null in step 7); none for pre / post
N_SHUFFLES <- if (MODEL %in% c("deathadj", "lead")) as.integer(Sys.getenv("N_SHUFFLES", "5")) else 0
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
LYMPHOID  <- c("Bcells", "CD4", "CD8", "NK")
DEATH     <- c("post_mortem", "ventilator")
OTHER_COV <- c("SMTSISCH", "SMRIN", "centre_B1", "male", "age") #true hk 11/14/2026

## ================= samples, proportions, covariates =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
# pre / post: keep only the pre- or post-mortem draws (replication check, summarised in step 7)
if (MODEL == "pre") cov_df <- cov_df[cov_df$post_mortem <= 0.5, ]
if (MODEL == "post") cov_df <- cov_df[cov_df$post_mortem > 0.5, ]
sample_ids <- cov_df$SAMPID
length(sample_ids)
props_file <- if (MODEL == "rawprops") "05_proportions_6classes_raw.csv" else "05_proportions_6classes.csv"
W_df <- read.csv(file.path(out_dir, props_file))
rownames(W_df) <- W_df$SAMPID
W <- as.matrix(W_df[sample_ids, CELLTYPES])

pcs <- fread(pc_file)
setnames(pcs, 1, "ID")
pcs <- pcs[ID %in% paste0("PC", 1:5)] #weird- i think its genotype
PC <- t(as.matrix(pcs[, cov_df$SUBJID, with = FALSE]))
colnames(PC) <- pcs$ID
varies <- function(v) {
  counts <- table(round(v, 6))
  length(counts) > 2 || (length(counts) == 2 && min(counts) >= 10)
}
death_used <- DEATH[sapply(DEATH, function(d) varies(cov_df[[d]]))]
death_used
death <- scale(as.matrix(cov_df[, death_used, drop = FALSE]))   # per cell type (C1)
C2 <- scale(cbind(as.matrix(cov_df[, OTHER_COV]), PC))         # bulk level (C2)
# 11.12: N_EPC = k (set above; main model 20): the top k hidden-noise PCs from step 5c (PCs of the expression left after
# removing composition + known covariates) as extra bulk covariates
if (N_EPC > 0) {
  epc <- read.csv(file.path(out_dir, "05c_expression_pcs.csv"))
  epc <- as.matrix(epc[match(sample_ids, epc$SAMPID), paste0("ePC", 1:N_EPC)])
  stopifnot(!anyNA(epc))
  C2 <- scale(cbind(C2, epc))
  rownames(C2) <- sample_ids
}
# 11.11: LABBATCH = 1 (robustness test, off by default): the GTEx expression batch (SMGEBTCH; it explains ~8% of the
# monocyte and B-cell proportions) as an extra bulk covariate. 132 batches of ~4 samples -> consecutive batches
# (LCSET number order = time) merged into groups of >= 30 samples, one indicator column per group (first group = baseline)
# Result 11.11 (chr21-22, OUT_TAG=_labbatch): hits nearly the same (Neut 52 -> 56, 49 in both; B 12 -> 13; CD8 8 -> 9;
# Mono 4 -> 4), effects at the same variant rho 0.94-1.00, ratio 0.96-1.01, all same sign -> not used in the main run.
# (The nucleic-acid batch SMNABTCH was not tested: 312 batches of ~2 samples, no proportion effect beyond permutation.)
                           #hk we dont adjust by lb batch i think in the end 11/14/2026

LABBATCH <- Sys.getenv("LABBATCH", "0") == "1"
if (LABBATCH) {
  ann <- fread(attr_ann_file, select = c("SAMPID", "SMGEBTCH"), quote = "")
  lab <- ann$SMGEBTCH[match(sample_ids, ann$SAMPID)]
  lcset_num <- as.numeric(sub("^LCSET-([0-9]+).*", "\\1", lab))          # "LCSET-10311,LCSET-10533" -> 10311
  stopifnot(!anyNA(lcset_num))
  batch_sizes <- table(lcset_num)                                           # samples per batch, in LCSET order
  group_of_batch <- integer(length(batch_sizes))
  g <- 1
  n_in_group <- 0
  for (k in seq_along(batch_sizes)) {                                       # fill a group until it has >= 30 samples
    group_of_batch[k] <- g
    n_in_group <- n_in_group + batch_sizes[k]
    if (n_in_group >= 30) {
      g <- g + 1
      n_in_group <- 0
    }
  }
  names(group_of_batch) <- names(batch_sizes)
  lab_group <- group_of_batch[as.character(lcset_num)]
  # merge a last group with < 30 samples into the one before
  last <- max(lab_group)
  if (sum(lab_group == last) < 30) lab_group[lab_group == last] <- last - 1
  cat("expression-batch groups:", length(unique(lab_group)), " samples per group:", range(table(lab_group)), "\n")
  lab_dummies <- model.matrix(~ factor(lab_group))[, -1, drop = FALSE]
  C2 <- scale(cbind(C2, lab_dummies))
}
rownames(death) <- sample_ids
rownames(C2) <- sample_ids
# W x each death variable (the per-cell-type death effects), as columns
W_death <- do.call(cbind, c(list(matrix(nrow = length(sample_ids), ncol = 0)),
                            lapply(seq_len(ncol(death)), function(k) W * death[, k])))

## ================= genes of this chromosome: expressed, unique symbol, TSS =================
gtf <- fread(cmd = paste("grep -P '\\tgene\\t'", gtf_file), header = FALSE)
gtf$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9)
gtf$tss <- ifelse(gtf$V7 == "+", gtf$V4, gtf$V5)
symbol_count <- table(gtf$gene_name)
unique_symbols <- names(symbol_count)[symbol_count == 1]
genes <- gtf %>% filter(V1 == CHR, gene_name %in% unique_symbols) %>% select(gene_name, tss)
genes <- genes %>% filter(gene_name %in% readLines(file.path(out_dir, "01_gene_subset.txt")))

tpm <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
X <- log2(tpm[genes$gene_name, sample_ids] + 1)
X <- X[apply(X, 1, sd) > 0, ]
genes <- genes %>% filter(gene_name %in% rownames(X))
if (N_TEST > 0) genes <- head(genes, N_TEST)
X <- X[genes$gene_name, , drop = FALSE]
# spike: only the genes of the step-8 plan
if (MODEL == "spike") {
  spike_plan <- read.csv(SPIKE_FILE)
  genes <- genes %>% filter(gene_name %in% spike_plan$gene_name)
  X <- X[genes$gene_name, , drop = FALSE]
}
nrow(genes)

## ================= a. TCA fit without genotype: sigma per gene and cell type, tau =================
null_file <- file.path(out_dir, paste0("06_scan_", CHR, model_tag, "_null_fit", if (N_TEST > 0) "_TEST", ".rds"))
# spike: reuse the main scan's null fit (same samples, same model; sigma per gene and cell type, tau)
if (MODEL %in% c("spike", "lead")) null_file <- file.path(out_dir, paste0("06_scan_", CHR, "_null_fit.rds"))
# LABBATCH: different covariates -> its own null fit
if (LABBATCH) null_file <- sub("_null_fit", "_labbatch_null_fit", null_file)
if (N_EPC > 0) null_file <- sub("_null_fit", paste0("_epc", N_EPC, "_null_fit"), null_file)
if (!file.exists(null_file)) {
  c1 <- if (ncol(death) > 0) death else NULL 
  m0 <- tca(X = X, W = W, C1 = c1, C2 = C2, parallel = N_CORES > 1, num_cores = N_CORES, verbose = FALSE) #this is actual tca :)
  saveRDS(list(sigmas = m0$sigmas_hat, tau = m0$tau), null_file)
}
null_fit <- readRDS(null_file)   # <- start here to skip
tau <- null_fit$tau
tau

## ================= genotypes: ALT dosage, MAF >= MIN_MAF =================
pl <- read.plink(paste0(geno_prefix, ".bed"), paste0(geno_prefix, ".bim"), paste0(geno_prefix, ".fam"))
G_all <- as(pl$genotypes, "numeric")                       # donors x variants: count of the bim's allele.2
alt <- sapply(strsplit(colnames(G_all), "_"), `[`, 4)
flip <- pl$map$allele.2 != alt
G_all[, flip] <- 2 - G_all[, flip]                         # now: count of the ALT allele
G_all <- G_all[match(cov_df$SUBJID, rownames(G_all)), ]
stopifnot(!anyNA(rownames(G_all)))
af <- colMeans(G_all, na.rm = TRUE) / 2
keep <- pmin(af, 1 - af) >= MIN_MAF
G_all <- G_all[, keep]
pos <- pl$map$position[keep]
missing <- which(is.na(G_all), arr.ind = TRUE)               # missing dosage -> the variant's mean
G_all[missing] <- colMeans(G_all, na.rm = TRUE)[missing[, 2]]
rownames(G_all) <- sample_ids
dim(G_all)

## ================= lead: the GTEx bulk lead variant of each gene =================
if (MODEL == "lead") {
  lead_tab <- fread(LEAD_FILE)
  lead_var <- setNames(lead_tab$variant_id, lead_tab$gene_name)
  lead_var <- lead_var[!duplicated(names(lead_var))]
  cat("genes with a GTEx lead variant that passed our MAF filter:",
      sum(lead_var[genes$gene_name] %in% colnames(G_all), na.rm = TRUE), "of", nrow(genes), "\n")
}

## ================= PHI = 1 (11.10, under test): genotype weights = RNA shares, not cell shares =================
# phi_c(i) = W_c(i) theta_c / sum_k W_k(i) theta_k (theta from step 5b). The genotype of cell type c enters as phi_c x g
# (on log2 expression, a change of log Z_c moves log bulk by phi_c x that change). Gate (as caNRD_editQTL): a cell type
# is tested in a gene only if its mean phi >= PHI_GATE; the others are left out of the model (their results = NA).
# The base design (W, W x death, covariates) is unchanged.
# 11.10: default ON (chr21-22 tests: minor-cell eGenes 0 -> 2-15 per cell type, calibration ~0.1 false per shuffle,
# 100% pre/post and raw-proportion replication, NK hits confirmed in OneK1K). PHI=0 = the old cell-share model.
PHI <- Sys.getenv("PHI", "1") == "1"
PHI_GATE <- 0.10 #true 
if (PHI) {
  theta_df <- read.csv(file.path(out_dir, "05b_theta.csv"))
  theta <- as.matrix(theta_df[, CELLTYPES])
  rownames(theta) <- theta_df$gene_name
  stopifnot(all(genes$gene_name %in% rownames(theta)))
}
phi_of <- function(gene) {                                   # samples x 6 RNA shares of one gene
  num <- sweep(W, 2, theta[gene, ], "*")
  num / rowSums(num)
}

## ================= spike: add a known effect in ONE cell type =================
# per planned gene: a random cis variant g (centred); in cell type c the expression rises by beta_sd x SD(gene) per
# allele, so the bulk (log2 expression) gets W_c(i) x g_i x beta_sd x SD(gene). Everything else stays real.
if (MODEL == "spike") {
  set.seed(2)
  spike_plan$variant_id <- NA_character_
  for (r in seq_len(nrow(spike_plan))) {
    j <- match(spike_plan$gene_name[r], genes$gene_name)
    cis <- which(abs(pos - genes$tss[j]) <= WINDOW)
    v <- cis[sample.int(length(cis), 1)]
    g <- G_all[, v] - mean(G_all[, v])
    ct <- spike_plan$celltype[r]
    # PHI = 1: the effect reaches the bulk through the cell type's RNA share (phi), the realistic way
    wt <- if (PHI) phi_of(genes$gene_name[j])[, ct] else W[, ct]
    X[genes$gene_name[j], ] <- X[genes$gene_name[j], ] + wt * g * spike_plan$beta_sd[r] * sd(X[genes$gene_name[j], ])
    spike_plan$variant_id[r] <- colnames(G_all)[v]
  }
  write.csv(spike_plan, file.path(out_dir, paste0("08_spike_plan_", CHR, "_with_variants", OUT_TAG, ".csv")), row.names = FALSE)
}
set.seed(1)
# N_SHUFFLES shuffles of the samples (genotype rows), each a full null scan; LD between variants is kept
shuffle_orders <- lapply(seq_len(N_SHUFFLES), function(k) sample(nrow(G_all)))
# pre / post: the main scan's best variant per gene and test (real genotype)
if (!MODEL %in% c("deathadj", "spike", "lead")) {
  # (11.10: the main scan with the same OUT_TAG, e.g. the PHI = 1 scan for PHI = 1 replication)
  main_variants <- read.csv(file.path(out_dir, paste0("06_scan_", CHR, "_top_variants", OUT_TAG, if (N_TEST > 0) "_TEST", ".csv")))
  main_variants <- main_variants[main_variants$set == "genotype", c("gene_name", "variant_id")]
}

## ================= b. one gene: weighted regression for every cis variant =================
lymph_idx <- match(LYMPHOID, CELLTYPES)
neut_idx <- match("Neutrophils", CELLTYPES)
contrast_nl <- numeric(6)
contrast_nl[neut_idx] <- 1
contrast_nl[lymph_idx] <- -1 / length(LYMPHOID)            # neutrophil minus mean lymphoid
het_C <- cbind(diag(5), -1)                                # effects 1..5 minus effect 6 (all equal <=> 0)

## ================= the regression for many variants (compiled kernel) =================
# 06_fit_variants.cpp does, per variant: genotype columns A[, a] * g, remove the base design, solve, robust HC3
# covariance, and returns the 6 effects, their SEs, the heterogeneity statistic and the neutrophil-lymphoid contrast.
# (Same math as the earlier R loop, results identical to ~1e-15; ~10x faster.) Here: those numbers -> p-values.
Rcpp::sourceCpp(file.path(pipe_dir, "06_fit_variants.cpp"))
fit_variants <- function(Gb, A, Q, y_r, h_base, df) {
  k <- fit_variants_cpp(Gb, A, Q, y_r, h_base, het_C, contrast_nl)
  b6 <- k[, 1:6, drop = FALSE]
  se <- k[, 7:12, drop = FALSE]
  p_ct <- 2 * pt(-abs(b6 / se), df)
  p_het <- pchisq(k[, 13], 5, lower.tail = FALSE)
  nl_est <- k[, 14]
  p_nl <- 2 * pt(-abs(nl_est / k[, 15]), df)
  out <- cbind(b6, se, p_ct, p_het, nl_est, p_nl)
  colnames(out) <- c(paste0("b_", CELLTYPES), paste0("se_", CELLTYPES), paste0("p_", CELLTYPES), "p_het", "nl_est", "p_nl")
  out
}
# PHI = 1: only the gated-in cell types (active, logical of length 6); results put back in the 6-column layout, NA for
# the others. het = are the active effects equal (k - 1 df); nl = neutrophil minus mean of the ACTIVE lymphoid types
# (NA if neutrophils or all lymphoid types are gated out).
fit_variants_phi <- function(Gb, A, Q, y_r, h_base, df, active) {
  k <- sum(active)
  hc <- if (k > 1) cbind(diag(k - 1), -1) else matrix(0, 0, k)
  cn <- numeric(k)
  act_ct <- CELLTYPES[active]
  if ("Neutrophils" %in% act_ct && any(LYMPHOID %in% act_ct)) {
    cn[act_ct == "Neutrophils"] <- 1
    cn[act_ct %in% LYMPHOID] <- -1 / sum(act_ct %in% LYMPHOID)
  }
  r <- fit_variants_k_cpp(Gb, A, Q, y_r, h_base, hc, cn, k)
  nv <- nrow(r)
  b6 <- matrix(NA_real_, nv, 6); se <- matrix(NA_real_, nv, 6)
  b6[, active] <- r[, seq_len(k)]
  se[, active] <- r[, k + seq_len(k)]
  p_ct <- 2 * pt(-abs(b6 / se), df)
  p_het <- if (k > 1) pchisq(r[, 2 * k + 1], k - 1, lower.tail = FALSE) else rep(NA_real_, nv)
  nl_est <- r[, 2 * k + 2]
  p_nl <- 2 * pt(-abs(nl_est / r[, 2 * k + 3]), df)
  out <- cbind(b6, se, p_ct, p_het, nl_est, p_nl)
  colnames(out) <- c(paste0("b_", CELLTYPES), paste0("se_", CELLTYPES), paste0("p_", CELLTYPES), "p_het", "nl_est", "p_nl")
  out
}

scan_gene <- function(j, shuffle = 0) {                     # shuffle = 0: real genotype; k: the k-th shuffle
  gene <- genes$gene_name[j]
  cis <- which(abs(pos - genes$tss[j]) <= WINDOW)
  # pre / post: test only the gene's best variants from the main (deathadj) scan, for replication
  if (!MODEL %in% c("deathadj", "spike", "lead")) cis <- cis[colnames(G_all)[cis] %in% main_variants$variant_id[main_variants$gene_name == gene]]
  if (MODEL == "lead") cis <- cis[colnames(G_all)[cis] %in% lead_var[gene]]
  if (length(cis) == 0) return(NULL)
  sig <- null_fit$sigmas[gene, ]
  sd_i <- sqrt(drop(W^2 %*% sig^2) + tau^2)                # expected SD of each sample
  y <- X[gene, ] / sd_i
  # base design, scaled by the sample SDs: W | W x each death variable | covariates
  base <- cbind(W, W_death, C2) / sd_i
  q <- qr(base)
  y_r <- qr.resid(q, y)
  df <- length(y) - q$rank - 6 - ncol(death)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  h_base <- rowSums(Q^2)                                    # leverage of each sample from the base design
  G <- G_all[, cis, drop = FALSE]
  if (shuffle > 0) G <- G[shuffle_orders[[shuffle]], , drop = FALSE]
  G <- sweep(G, 2, colMeans(G))                            # centred genotypes, samples x variants
  A <- cbind(W, death) / sd_i                              # the genotype enters as A[, a] * g (a = 1..6 cell types, then death)
  # all variants of the gene at once, in blocks of 1000 (memory)
  blocks <- split(seq_len(ncol(G)), ceiling(seq_len(ncol(G)) / 1000))
  if (!PHI) {
    res <- do.call(rbind, lapply(blocks, function(idx) {
      fit_variants(G[, idx, drop = FALSE], A, Q, y_r, h_base, df)
    }))
  } else {
    # PHI = 1: genotype columns = RNA shares of the gated-in cell types (then death)
    phi <- phi_of(gene)
    active <- colMeans(phi) >= PHI_GATE
    if (!any(active)) return(NULL)
    A <- cbind(phi[, active, drop = FALSE], death) / sd_i
    df <- length(y) - q$rank - sum(active) - ncol(death)
    res <- do.call(rbind, lapply(blocks, function(idx) {
      fit_variants_phi(G[, idx, drop = FALSE], A, Q, y_r, h_base, df, active)
    }))
  }
  res <- data.frame(gene_name = gene, variant_id = colnames(G_all)[cis], res)
  # keep variants with at least one result (PHI: gated cell types and het can be NA)
  res <- if (PHI) res[rowSums(!is.na(res[, paste0("p_", CELLTYPES)])) > 0, ] else res[!is.na(res$p_het), ]
  if (nrow(res) == 0) return(NULL)
  # keep each test's best variant, and the number of variants (for the Bonferroni)
  tests <- c(paste0("p_", CELLTYPES), "p_het", "p_nl")
  tests <- tests[sapply(tests, function(tt) any(!is.na(res[[tt]])))]   # PHI: skip tests that are all NA
  top <- bind_rows(lapply(tests, function(tt) {
    best <- res[which.min(res[[tt]]), ]
    best$test <- sub("^p_", "", tt)
    best$p_min <- best[[tt]]
    best$n_variants <- nrow(res)
    best
  }))
  top$set <- if (shuffle > 0) paste0("shuffled_", shuffle) else "genotype"
  top
}

## ================= run: all genes, real and shuffled =================
idx <- seq_len(nrow(genes))
# 11.09: a parallel worker killed for memory returns NULL or an error, which bind_rows silently dropped
# (genome-wide run lost 1-2 genes on chr1/2/4). So: run in parallel, then refit every gene that came back
# empty or failed one at a time (genes with no testable variant just return NULL again, quickly),
# and stop if any gene still fails.
scan_all_genes <- function(shuffle) {
  res <- mclapply(idx, scan_gene, shuffle = shuffle, mc.cores = N_CORES, mc.preschedule = FALSE)
  # genes to refit: no result, or an error
  redo <- which(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1)))
  res[redo] <- lapply(redo, function(j) try(scan_gene(j, shuffle = shuffle), silent = TRUE))
  # genes still failing: stop
  failed <- which(vapply(res, function(r) inherits(r, "try-error"), logical(1)))
  if (length(failed) > 0) stop("scan_gene failed for ", length(failed), " genes (shuffle ", shuffle, ")")
  cat("shuffle", shuffle, ": refitted", length(redo), "empty genes one by one;",
      sum(vapply(res, is.null, logical(1))), "genes have no testable variant\n")
  bind_rows(res)
}
top <- scan_all_genes(0)
for (k in seq_len(N_SHUFFLES)) {
  shuffled <- scan_all_genes(k)
  top <- bind_rows(top, shuffled)
}
top$chr <- CHR

## ================= gene level: Bonferroni within gene, BH across genes =================
# (theoretical p-values; the calibrated empirical p-values against the pooled shuffles are made in step 7)
top$p_gene <- pmin(top$p_min * top$n_variants, 1)
top <- top %>% group_by(set, test) %>% mutate(fdr = p.adjust(p_gene, "BH")) %>% ungroup()
summary_tab <- top %>% group_by(test, set) %>% summarise(genes = n(), egenes_fdr05 = sum(fdr < 0.05), .groups = "drop")
print(as.data.frame(summary_tab))

## ================= save =================
top$model <- MODEL
# OUT_TAG (11.10): extra suffix for test runs (e.g. "_shuf25" = more shuffles), so the main outputs are not overwritten
suffix <- paste0(model_tag, OUT_TAG, if (N_TEST > 0) "_TEST" else "")
gene_cols <- c("chr", "model", "gene_name", "test", "set", "variant_id", "n_variants", "p_min", "p_gene", "fdr")
write.csv(top[, gene_cols], file.path(out_dir, paste0("06_scan_", CHR, "_genes", suffix, ".csv")), row.names = FALSE)
write.csv(top, file.path(out_dir, paste0("06_scan_", CHR, "_top_variants", suffix, ".csv")), row.names = FALSE)
message("06_cis_scan done: ", CHR, " ", MODEL)
#verified on 11/14/2026 by hk
