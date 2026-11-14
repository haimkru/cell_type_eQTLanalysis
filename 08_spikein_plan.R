#author: haim krupkin+claude
#date: 11.10.2026
#description: step 8a. Power check (spike-in) plan. Genome-wide, only neutrophils had calibrated eGenes. Is that because
# the minor cell types have no power at n = 670, or because the pipeline cannot find them? To answer, we ADD known
# cell-type eQTLs to real expression and see whether step 6 finds them (step 8b scores it).
# Plan: genes of one chromosome (CHR, default chr22) with NO signal in the main scan (smallest gene-level p of all
# tests > 0.1). Each gets one cell type and one effect size (beta_sd = effect per allele in that cell type, in SDs of
# the gene's log2 expression). 10 genes per cell type x effect size, plus 20 genes with beta_sd = 0 (null control).
# Step 6 (MODEL=spike) picks a random cis variant per gene and adds W_c x genotype x beta_sd x SD(gene).
# Input : output/06_scan_<CHR>_genes.csv (main scan)
# Output: output/08_spike_plan_<CHR>.csv (gene_name, celltype, beta_sd)
suppressMessages(library(dplyr))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHR      <- Sys.getenv("CHR", "chr22")
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
BETAS     <- c(0.5, 1, 2, 4)        # per-allele effect in the cell type, in SD of the gene
N_PER     <- 10                     # genes per cell type x effect size
N_NULL    <- 20                     # genes with no added effect (beta_sd = 0)

# 11.10: SCAN_TAG = which main scan defines "no real signal" (e.g. "_phi_win100k_shuf25"); PHI = 1: a gene is planted
# only in a cell type with mean RNA share (phi, step 5b) >= MIN_PHI (a cell type cannot carry an eQTL in a gene it
# does not express)
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")   # 11.13: main result = 20 hidden-noise PCs; "" = the model without
PHI      <- Sys.getenv("PHI", "1") == "1"                  # default ON, as step 6
MIN_PHI  <- 0.2

## ================= genes with no real signal =================
main <- read.csv(file.path(out_dir, paste0("06_scan_", CHR, "_genes", SCAN_TAG, ".csv")))
quiet <- main %>%
  filter(set == "genotype") %>%
  group_by(gene_name) %>%
  summarise(p_best = min(p_gene), .groups = "drop") %>%
  filter(p_best > 0.1)
nrow(quiet)
if (PHI) {
  mean_phi <- read.csv(file.path(out_dir, "05b_mean_phi.csv"))
  quiet <- quiet %>% left_join(mean_phi, by = "gene_name")
}

## ================= the plan =================
# genes drawn without replacement; with PHI, each cell type draws from the quiet genes it expresses (phi >= MIN_PHI)
set.seed(1)
plan <- expand.grid(rep = seq_len(N_PER), beta_sd = BETAS, celltype = CELLTYPES, stringsAsFactors = FALSE)
plan <- bind_rows(plan, data.frame(rep = seq_len(N_NULL), beta_sd = 0, celltype = "Neutrophils"))
plan$gene_name <- NA_character_
used <- character(0)
for (ct in c(CELLTYPES[CELLTYPES != "Neutrophils"], "Neutrophils")) {   # rarest-expressed first
  rows <- which(plan$celltype == ct)
  pool <- setdiff(quiet$gene_name[if (PHI) quiet[[ct]] >= MIN_PHI else TRUE], used)
  stopifnot(length(pool) >= length(rows))
  plan$gene_name[rows] <- sample(pool, length(rows))
  used <- c(used, plan$gene_name[rows])
}
plan <- plan %>% select(gene_name, celltype, beta_sd)
table(plan$celltype, plan$beta_sd)
write.csv(plan, file.path(out_dir, paste0("08_spike_plan_", CHR, SCAN_TAG, ".csv")), row.names = FALSE)
