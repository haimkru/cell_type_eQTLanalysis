#author: haim krupkin+claude
#date: 11.10.2026
#description: step 5b. Expression of each gene in each cell type (theta), for the RNA-share weights (phi) of step 6
# (MODEL=phi, under test). Why: step 6 weighted each cell type's genotype effect by its share of CELLS (W). A cell type's
# effect on bulk expression is really weighted by its share of the gene's RNA:
#   phi_c(i) = W_c(i) * theta_c / sum_k W_k(i) * theta_k
# e.g. a B-cell gene: B cells are ~5% of cells but ~100% of its RNA. caNRD_editQTL uses the same phi for edQTLs.
# theta per gene: non-negative least squares of the linear TPM on the 6 proportions over all samples
# (TPM_i = sum_c W_c(i) theta_c), floor 1e-3 (= "not expressed"), as caNRD_editQTL's GTEx theta.
# Inputs : output/01_tpm_noHb.rds (linear TPM, Hb-free), 01_gene_subset.txt, 05_proportions_6classes.csv, 02_covariates.csv
# Output : output/05b_theta.csv (gene_name + 6 cell types), output/05b_mean_phi.csv (mean phi per gene and cell type)
suppressMessages(library(dplyr))
suppressMessages(library(nnls))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
THETA_FLOOR <- 1e-3

## ================= inputs =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
sample_ids <- cov_df$SAMPID
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
rownames(W_df) <- W_df$SAMPID
W <- as.matrix(W_df[sample_ids, CELLTYPES])
genes <- readLines(file.path(out_dir, "01_gene_subset.txt"))
tpm <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
genes <- intersect(genes, rownames(tpm))
length(genes)
Y <- tpm[genes, sample_ids]

## ================= theta: NNLS per gene =================
theta <- t(apply(Y, 1, function(y) nnls(W, as.numeric(y))$x))
colnames(theta) <- CELLTYPES
theta <- pmax(theta, THETA_FLOOR)
cat("share of theta entries at the floor:", round(mean(theta == THETA_FLOOR), 3), "\n")

## ================= mean phi per gene and cell type =================
mean_phi <- t(sapply(seq_len(nrow(theta)), function(j) {
  num <- sweep(W, 2, theta[j, ], "*")
  colMeans(num / rowSums(num))
}))
colnames(mean_phi) <- CELLTYPES
cat("genes with mean phi >= 0.10, per cell type:\n")
print(colSums(mean_phi >= 0.10))

write.csv(data.frame(gene_name = genes, theta, check.names = FALSE), file.path(out_dir, "05b_theta.csv"), row.names = FALSE)
write.csv(data.frame(gene_name = genes, mean_phi, check.names = FALSE), file.path(out_dir, "05b_mean_phi.csv"), row.names = FALSE)
message("05b_theta done")
