#author: haim krupkin+claude
#date: 11.14.2026
#description: step 12. Cell-type TRANS-eQTL scan: every candidate variant (step 12a) x every expressed gene of one
# chromosome (CHR) that is FAR from the variant (other chromosome, or > TRANS_DIST on the same one).
# The model is exactly the main cis model of step 6 (deathadj, phi weights + phi gate 0.10, 20 hidden-noise PCs, robust
# HC3 errors, same compiled kernel): per gene, the 6 cell-type genotype effects (gamma_h) are estimated together, plus
# the heterogeneity test (het: are they equal?) and neutrophil vs mean lymphoid (nl). Only the variant set changes.
# The per-gene noise model (sigma_h, tau) is the one step 6 already fitted for that gene (06_scan_<CHR>_epc20_null_fit.rds),
# so this is the cis engine pointed at distant variants, nothing refitted.
# Phi gate: cell type h is tested for a gene only if it makes >= 10% of that gene's RNA (as in cis): a trans effect
# "in cell type h" is a change of THIS gene's expression in h.
# Null: the same scan with the genotype rows shuffled (N_SHUFFLES, same shuffle seed as step 6). Shuffling whole rows
# keeps the LD between candidates and the correlation between genes, so the null has the same dependence as the real scan.
# Storage: ~2,700 variants x ~19,000 genes x 8 tests per set is too much to keep, so per gene and set only
#   - pairs with any test p < P_KEEP (all effects, SEs, p-values), and
#   - the number of tests done per test type (for the FDR denominator)
# are saved. Calls (step 12b) need p far below P_KEEP (genome-wide trans: ~1e-9), so nothing needed is lost.
# Inputs : output/12a_trans_candidates<SCAN_TAG>.csv, output/genotypes/<chr>.{bed,bim,fam}, the step-6 inputs
# Outputs: output/12_trans_<CHR>_pairs<OUT_TAG>.csv   (kept pairs, real + shuffled; column set)
#          output/12_trans_<CHR>_ntests<OUT_TAG>.csv  (set x test: number of variant-gene tests done)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))
library(parallel)

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHR      <- Sys.getenv("CHR", "chr22")                     # chromosome of the GENES (variants: all chromosomes)
gtex     <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8"
gtf_file <- file.path(gtex, "references/gencode.v26.GRCh38.genes.gtf")
pc_file  <- file.path(gtex, "eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt")
TRANS_DIST <- as.numeric(Sys.getenv("TRANS_DIST", "5e6"))  # same chromosome: variant > 5 Mb from the TSS = trans
P_KEEP   <- as.numeric(Sys.getenv("P_KEEP", "1e-4"))      # keep pairs with any p below this
N_CORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
N_TEST   <- as.integer(Sys.getenv("N_TEST", "0"))          # > 0: only the first N_TEST genes (quick test)
N_SHUFFLES <- as.integer(Sys.getenv("N_SHUFFLES", "5"))
N_EPC    <- as.integer(Sys.getenv("N_EPC", "20"))          # must match the step-6 run whose null fit is reused
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")               # tag of the step-6 run (null fit) and of step 12a
CAND_FILE <- Sys.getenv("CAND_FILE", file.path(out_dir, paste0("12a_trans_candidates", SCAN_TAG, ".csv")))
OUT_TAG  <- Sys.getenv("OUT_TAG", SCAN_TAG)
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
LYMPHOID  <- c("Bcells", "CD4", "CD8", "NK")
DEATH     <- c("post_mortem", "ventilator")
OTHER_COV <- c("SMTSISCH", "SMRIN", "centre_B1", "male", "age")
PHI_GATE  <- 0.10

## ================= samples, proportions, covariates (as step 6, MODEL = deathadj) =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
sample_ids <- cov_df$SAMPID
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
rownames(W_df) <- W_df$SAMPID
W <- as.matrix(W_df[sample_ids, CELLTYPES])
# genotype PCs 1-5 (GTEx)
pcs <- fread(pc_file)
setnames(pcs, 1, "ID")
pcs <- pcs[ID %in% paste0("PC", 1:5)]
PC <- t(as.matrix(pcs[, cov_df$SUBJID, with = FALSE]))
colnames(PC) <- pcs$ID
# death variables per cell type (C1), other covariates + PCs at bulk level (C2), hidden-noise PCs
death <- scale(as.matrix(cov_df[, DEATH, drop = FALSE]))
C2 <- scale(cbind(as.matrix(cov_df[, OTHER_COV]), PC))
if (N_EPC > 0) {
  epc <- read.csv(file.path(out_dir, "05c_expression_pcs.csv"))
  epc <- as.matrix(epc[match(sample_ids, epc$SAMPID), paste0("ePC", 1:N_EPC)])
  stopifnot(!anyNA(epc))
  C2 <- scale(cbind(C2, epc))
}
rownames(death) <- sample_ids
rownames(C2) <- sample_ids
# W x each death variable (per-cell-type death effects)
W_death <- do.call(cbind, lapply(seq_len(ncol(death)), function(k) W * death[, k]))

## ================= genes of CHR: the genes step 6 fitted (its null fit), with TSS =================
null_fit <- readRDS(file.path(out_dir, paste0("06_scan_", CHR, SCAN_TAG, "_null_fit.rds")))   # <- start here to skip
tau <- null_fit$tau
gtf <- fread(cmd = paste("grep -P '\\tgene\\t'", gtf_file), header = FALSE)
gtf$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9)
gtf$tss <- ifelse(gtf$V7 == "+", gtf$V4, gtf$V5)
genes <- gtf %>%
  filter(V1 == CHR, gene_name %in% rownames(null_fit$sigmas)) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  select(gene_name, tss)
if (N_TEST > 0) genes <- head(genes, N_TEST)
tpm <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
X <- log2(tpm[genes$gene_name, sample_ids, drop = FALSE] + 1)   # log2(TPM + 1), Hb-free (step 1)
rm(tpm)
cat(CHR, ": genes", nrow(genes), "\n")

## ================= phi (RNA share per cell type), as step 6 =================
theta_df <- read.csv(file.path(out_dir, "05b_theta.csv"))
theta <- as.matrix(theta_df[, CELLTYPES])
rownames(theta) <- theta_df$gene_name
stopifnot(all(genes$gene_name %in% rownames(theta)))
phi_of <- function(gene) {                                   # samples x 6 RNA shares of one gene
  num <- sweep(W, 2, theta[gene, ], "*")
  num / rowSums(num)
}

## ================= candidate genotypes: ALT dosage, all chromosomes =================
cand <- read.csv(CAND_FILE)
read_candidates <- function(ch) {
  ids <- cand$variant_id[cand$chr == ch]
  pfx <- file.path(out_dir, "genotypes", ch)
  pl <- read.plink(paste0(pfx, ".bed"), paste0(pfx, ".bim"), paste0(pfx, ".fam"), select.snps = ids)
  G <- as(pl$genotypes, "numeric")                           # count of the bim's allele.2
  alt <- sapply(strsplit(colnames(G), "_"), `[`, 4)
  flip <- pl$map$allele.2 != alt
  G[, flip] <- 2 - G[, flip]                                 # count of the ALT allele
  G[match(cov_df$SUBJID, rownames(G)), , drop = FALSE]
}
G_all <- do.call(cbind, lapply(unique(cand$chr), read_candidates))
stopifnot(!anyNA(rownames(G_all)))
missing <- which(is.na(G_all), arr.ind = TRUE)               # missing dosage -> the variant's mean
G_all[missing] <- colMeans(G_all, na.rm = TRUE)[missing[, 2]]
rownames(G_all) <- sample_ids
G_all <- sweep(G_all, 2, colMeans(G_all))                    # centred
var_chr <- sub("_.*", "", colnames(G_all))
var_pos <- as.numeric(sapply(strsplit(colnames(G_all), "_"), `[`, 2))
cat("candidate variants:", ncol(G_all), "\n")
set.seed(1)                                                  # the same shuffle orders as step 6
shuffle_orders <- lapply(seq_len(N_SHUFFLES), function(k) sample(nrow(G_all)))

## ================= the per-variant regression (compiled kernel of step 6) =================
Rcpp::sourceCpp(file.path(pipe_dir, "06_fit_variants.cpp"))
TESTS <- c(paste0("p_", CELLTYPES), "p_het", "p_nl")
fit_block <- function(Gb, A, Q, y_r, h_base, df, active) {
  k <- sum(active)
  hc <- if (k > 1) cbind(diag(k - 1), -1) else matrix(0, 0, k)
  cn <- numeric(k)
  act_ct <- CELLTYPES[active]
  # neutrophil minus the mean of the ACTIVE lymphoid types (none -> no contrast)
  if ("Neutrophils" %in% act_ct && any(LYMPHOID %in% act_ct)) {
    cn[act_ct == "Neutrophils"] <- 1
    cn[act_ct %in% LYMPHOID] <- -1 / sum(act_ct %in% LYMPHOID)
  }
  r <- fit_variants_k_cpp(Gb, A, Q, y_r, h_base, hc, cn, k)
  nv <- nrow(r)
  b6 <- matrix(NA_real_, nv, 6)
  se <- matrix(NA_real_, nv, 6)
  b6[, active] <- r[, seq_len(k)]
  se[, active] <- r[, k + seq_len(k)]
  p_ct <- 2 * pt(-abs(b6 / se), df)
  p_het <- if (k > 1) pchisq(r[, 2 * k + 1], k - 1, lower.tail = FALSE) else rep(NA_real_, nv)
  nl_est <- r[, 2 * k + 2]
  p_nl <- 2 * pt(-abs(nl_est / r[, 2 * k + 3]), df)
  out <- cbind(b6, se, p_ct, p_het, nl_est, p_nl)
  colnames(out) <- c(paste0("b_", CELLTYPES), paste0("se_", CELLTYPES), TESTS[1:6], "p_het", "nl_est", "p_nl")
  out
}

## ================= one gene: all its trans candidates, real and shuffled =================
scan_gene <- function(j) {
  gene <- genes$gene_name[j]
  # trans = other chromosome, or same chromosome and > TRANS_DIST from the TSS
  trans <- which(var_chr != CHR | abs(var_pos - genes$tss[j]) > TRANS_DIST)
  # phi gate: the cell types making >= 10% of this gene's RNA
  phi <- phi_of(gene)
  active <- colMeans(phi) >= PHI_GATE
  if (!any(active) || length(trans) == 0) return(NULL)
  # weighted design, as step 6: each sample scaled by its expected SD
  sig <- null_fit$sigmas[gene, ]
  sd_i <- sqrt(drop(W^2 %*% sig^2) + tau^2)
  y <- X[gene, ] / sd_i
  base <- cbind(W, W_death, C2) / sd_i
  q <- qr(base)
  y_r <- qr.resid(q, y)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  h_base <- rowSums(Q^2)
  A <- cbind(phi[, active, drop = FALSE], death) / sd_i
  df <- length(y) - q$rank - sum(active) - ncol(death)
  blocks <- split(trans, ceiling(seq_along(trans) / 1000))
  pairs <- list()
  ntests <- list()
  for (s in 0:N_SHUFFLES) {
    rows <- if (s == 0) seq_len(nrow(G_all)) else shuffle_orders[[s]]
    res <- do.call(rbind, lapply(blocks, function(idx) {
      fit_block(G_all[rows, idx, drop = FALSE], A, Q, y_r, h_base, df, active)
    }))
    set_name <- if (s == 0) "genotype" else paste0("shuffled_", s)
    # number of tests done per test type (non-NA p)
    ntests[[set_name]] <- data.frame(set = set_name, test = sub("^p_", "", TESTS),
                                     n = colSums(!is.na(res[, TESTS, drop = FALSE])))
    # keep pairs with any p < P_KEEP
    keep <- which(apply(res[, TESTS, drop = FALSE], 1, function(p) any(p < P_KEEP, na.rm = TRUE)))
    if (length(keep) > 0) {
      pairs[[set_name]] <- data.frame(gene_name = gene, gene_chr = CHR, gene_tss = genes$tss[j],
                                      variant_id = colnames(G_all)[trans[keep]], set = set_name,
                                      res[keep, , drop = FALSE])
    }
  }
  list(pairs = bind_rows(pairs), ntests = bind_rows(ntests))
}

## ================= run all genes (failed workers refitted one by one, as step 6) =================
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(genes)), function(j) try(scan_gene(j), silent = TRUE),
                mc.cores = N_CORES, mc.preschedule = FALSE)
redo <- which(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1)))
res[redo] <- lapply(redo, function(j) try(scan_gene(j), silent = TRUE))
failed <- which(vapply(res, function(r) inherits(r, "try-error"), logical(1)))
if (length(failed) > 0) stop("scan_gene failed for ", length(failed), " genes, e.g. ", res[[failed[1]]])
cat("genes refitted one by one:", length(redo), "; time:", format(Sys.time() - t0), "\n")
res <- res[!vapply(res, is.null, logical(1))]
pairs <- bind_rows(lapply(res, `[[`, "pairs"))
ntests <- bind_rows(lapply(res, `[[`, "ntests")) %>%
  group_by(set, test) %>%
  summarise(n_tests = sum(n), .groups = "drop") %>%
  mutate(chr = CHR)
print(as.data.frame(ntests))
# kept pairs per set (real should exceed shuffled if there is signal)
print(table(pairs$set))

## ================= save =================
suffix <- paste0(OUT_TAG, if (N_TEST > 0) "_TEST" else "")
pairs_file <- file.path(out_dir, paste0("12_trans_", CHR, "_pairs", suffix, ".csv"))
ntests_file <- file.path(out_dir, paste0("12_trans_", CHR, "_ntests", suffix, ".csv"))
fwrite(pairs, pairs_file)
fwrite(ntests, ntests_file)
c(pairs_file, ntests_file)
message("12_trans_scan done: ", CHR)
