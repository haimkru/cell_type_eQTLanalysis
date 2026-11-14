#author: haim krupkin+claude
#date: 11.12.2026
#description: step 5c (under test, 11.12). Hidden-noise factors for the eQTL model: principal components of the
# expression that is LEFT after removing cell composition and the known covariates.
# Why not plain expression PCs (or GTEx PEER factors)? In blood those mostly capture cell composition (neutrophil vs
# lymphocyte ratio), which is exactly the variation step 6 learns the cell-type effects from. Removing composition
# first keeps only hidden noise (library prep, RNA quality beyond RIN, unknown batches).
# How:
#   1. log2(TPM + 1) of the expressed genes (step 1 subset), all 670 samples
#   2. per gene: regress on the 6 proportions, proportions x death (blood draw, ventilator), ischemic time, RIN,
#      centre, sex, age, 5 genotype PCs (= the step-6 design without genotype); keep the residuals
#   3. PCA of the residuals (each gene scaled to SD 1); keep the top 20 PCs
# Inputs : output/01_tpm_noHb.rds, 01_gene_subset.txt, 02_covariates.csv, 05_proportions_6classes.csv, GTEx genotype PCs
# Outputs: output/05c_expression_pcs.csv (SAMPID + ePC1..ePC20), output/05c_variance_explained.csv
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
pc_file  <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
DEATH     <- c("post_mortem", "ventilator")
OTHER_COV <- c("SMTSISCH", "SMRIN", "centre_B1", "male", "age")
N_PCS     <- 20

## ================= read =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
sample_ids <- cov_df$SAMPID
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
W <- as.matrix(W_df[match(sample_ids, W_df$SAMPID), CELLTYPES])
pcs <- fread(pc_file)
setnames(pcs, 1, "ID")
geno_pc <- t(as.matrix(pcs[ID %in% paste0("PC", 1:5), cov_df$SUBJID, with = FALSE]))

tpm <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
genes <- intersect(readLines(file.path(out_dir, "01_gene_subset.txt")), rownames(tpm))
Y <- t(log2(tpm[genes, sample_ids] + 1))                     # samples x genes
Y <- Y[, apply(Y, 2, sd) > 0]
dim(Y)

## ================= residuals after composition + known covariates =================
# W sums to 1 per sample, so it already contains the intercept (no separate intercept column)
death <- as.matrix(cov_df[, DEATH])
W_death <- do.call(cbind, lapply(seq_len(ncol(death)), function(k) W * death[, k]))
design <- cbind(W, W_death, scale(as.matrix(cov_df[, OTHER_COV])), scale(geno_pc))
fit <- lm.fit(design, Y)                                     # all genes at once
R <- fit$residuals
cat("design columns:", ncol(design), " rank:", fit$rank, "\n")

## ================= PCA of the residuals =================
p <- prcomp(R, center = TRUE, scale. = TRUE, rank. = N_PCS)
var_expl <- (p$sdev^2 / sum(p$sdev^2))[1:N_PCS]
print(round(var_expl, 3))
# check: the PCs should be ~uncorrelated with the proportions (composition was removed)
cat("max |cor| of an ePC with a proportion:", round(max(abs(cor(p$x, W))), 3), "\n")

ePC <- data.frame(SAMPID = sample_ids, p$x)
names(ePC) <- c("SAMPID", paste0("ePC", 1:N_PCS))
write.csv(ePC, file.path(out_dir, "05c_expression_pcs.csv"), row.names = FALSE)
write.csv(data.frame(pc = paste0("ePC", 1:N_PCS), variance_explained = var_expl),
          file.path(out_dir, "05c_variance_explained.csv"), row.names = FALSE)
message("05c_expression_pcs done")
