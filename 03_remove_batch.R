#author: haim krupkin+claude
#date: 11.07.2026
#description: step 3. Remove the way-of-death batch effect from the expression BEFORE CIBERSORTx.
# Why: blood draw (pre/post-mortem) explains ~70% of expression PC1 and shifts the estimated proportions
# (neutrophils 0.64 pre vs 0.36 post). We want proportions that vary between people, not between ways of death.
# How (linear scale, so CIBERSORTx's "bulk = sum of cell types" still holds):
#   per gene, a log-link regression on the TPM itself (quasi-Poisson GLM):
#       E[TPM] = exp(b0 + post_mortem + ventilator + ischemic time + RIN + centre B1)
#   then each value is divided by its fitted batch fold change:  TPM_new = TPM / exp(batch part)
#   -> never negative, zeros stay zeros, the gene's average level is kept.
# Sex and age are NOT removed. Genes with < 10 non-zero samples are left as they are. Then re-TPM.
# Inputs : output/01_tpm_noHb.rds, output/02_covariates.csv
# Outputs: output/03_tpm_noHb_batch_removed.tsv (CIBERSORTx mixture), output/03_batch_check.csv
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
library(parallel)

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
BATCH_VARS <- c("post_mortem", "ventilator", "SMTSISCH", "SMRIN", "centre_B1")
N_CORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

## ================= read =================
X <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
rownames(cov_df) <- cov_df$SAMPID
cov_df <- cov_df[colnames(X), ]
B <- scale(as.matrix(cov_df[, BATCH_VARS]))      # batch variables, centred and scaled
D <- cbind(intercept = 1, B)                     # design matrix

## ================= per gene: fit, then divide by the batch fold change =================
remove_batch <- function(x) {
  if (sum(x > 0) < 10) return(x)                                        # (almost) never expressed
  f <- tryCatch(glm.fit(D, x, family = quasipoisson()), error = function(e) NULL)
  if (is.null(f)) return(x)
  if (!f$converged) return(x)
  batch_part <- drop(B %*% f$coefficients[-1])                          # log fold change from the batch variables
  x / exp(batch_part)
}
rows <- mclapply(seq_len(nrow(X)), function(i) remove_batch(X[i, ]), mc.cores = N_CORES)
X_new <- do.call(rbind, rows)
dimnames(X_new) <- dimnames(X)

# re-TPM
X_new <- sweep(X_new, 2, colSums(X_new), "/")
X_new <- X_new * 1e6
stopifnot(min(X_new) >= 0)

saveRDS(X_new, file.path(out_dir, "03_tpm_noHb_batch_removed.rds"))
X_new <- readRDS(file.path(out_dir, "03_tpm_noHb_batch_removed.rds"))   # <- start here to skip

## ================= check: does blood draw still explain PC1? =================
pc1_adj_r2 <- function(M) {
  L <- log2(M + 1)
  gene_var <- apply(L, 1, var)
  top <- order(-gene_var)[1:5000]                                       # 5,000 most variable genes
  pc <- prcomp(t(L[top, ]), center = TRUE, scale. = TRUE)
  fit <- lm(pc$x[, 1] ~ cov_df$post_mortem)
  summary(fit)$adj.r.squared
}
check <- data.frame(data = c("before", "after"),
                    PC1_adj_r2_blood_draw = c(pc1_adj_r2(X), pc1_adj_r2(X_new)))
print(check)
write.csv(check, file.path(out_dir, "03_batch_check.csv"), row.names = FALSE)

## ================= save the CIBERSORTx mixture =================
tsv_out <- data.table(GeneSymbol = rownames(X_new), X_new)
out_file <- file.path(out_dir, "03_tpm_noHb_batch_removed.tsv")
fwrite(tsv_out, out_file, sep = "\t")
out_file
message("03_remove_batch done")
