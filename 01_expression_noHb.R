#author: haim krupkin+claude
#date: 11.07.2026
#description: step 1. GTEx whole-blood expression without haemoglobin genes.
#   1. read the GTEx TPM (670 whole-blood samples)
#   2. remove the 10 haemoglobin genes (they are > half of all reads in blood and distort every per-sample scaling)
#   3. re-TPM: each sample sums to 1e6 again
#   4. add up rows that share a gene symbol (CIBERSORTx works on symbols)
#   5. gene subset = genes with TPM > 0.1 in >= 20% of samples (used later for the eQTL step)
# Inputs : expression/gtex_wholeblood_gene_tpm.csv.gz
# Outputs: output/01_tpm_noHb_all_genes.tsv  (linear TPM, all genes; first column GeneSymbol)
#          output/01_gene_subset.txt         (expressed genes)
#          output/01_tpm_noHb.rds            (same matrix as an R object)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
tpm_file <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/expression/gtex_wholeblood_gene_tpm.csv.gz"
HB_GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "HBE1")
TPM_THRESHOLD <- 0.1      # gene subset: TPM above this ...
MIN_FRACTION  <- 0.2      # ... in at least this share of samples

## ================= 1. read =================
tpm_raw <- fread(tpm_file)
samples <- setdiff(colnames(tpm_raw), c("Name", "Description"))
length(samples)

## ================= 2. remove haemoglobin genes =================
is_hb <- tpm_raw$Description %in% HB_GENES
sum(is_hb)
tpm_noHb <- tpm_raw[!is_hb, ]

## ================= 3. re-TPM =================
tpm_mat <- as.matrix(tpm_noHb[, ..samples])
sample_totals <- colSums(tpm_mat)
tpm_mat <- sweep(tpm_mat, 2, sample_totals, "/")   # each sample now sums to 1
tpm_mat <- tpm_mat * 1e6                            # ... and now to 1e6
stopifnot(all(abs(colSums(tpm_mat) - 1e6) < 1e-3))

## ================= 4. sum rows with the same gene symbol =================
mix_all <- rowsum(tpm_mat, group = tpm_noHb$Description, reorder = FALSE)
dim(mix_all)
stopifnot(!any(HB_GENES %in% rownames(mix_all)))

## ================= 5. gene subset (expressed genes) =================
min_samples <- ceiling(MIN_FRACTION * length(samples))
n_above <- rowSums(mix_all > TPM_THRESHOLD)
gene_subset <- rownames(mix_all)[n_above >= min_samples]
length(gene_subset)

## ================= save =================
saveRDS(mix_all, file.path(out_dir, "01_tpm_noHb.rds"))
mix_all <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))   # <- start here to skip

tsv_out <- data.table(GeneSymbol = rownames(mix_all), mix_all)
fwrite(tsv_out, file.path(out_dir, "01_tpm_noHb_all_genes.tsv"), sep = "\t")
writeLines(gene_subset, file.path(out_dir, "01_gene_subset.txt"))

out_files <- file.path(out_dir, c("01_tpm_noHb_all_genes.tsv", "01_gene_subset.txt", "01_tpm_noHb.rds"))
out_files
stopifnot(all(file.exists(out_files)))
message("01_expression_noHb done")
