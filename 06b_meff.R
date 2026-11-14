#author: haim krupkin+claude
#date: 11.10.2026
#description: step 6b. Effective number of independent tests (M_eff) per gene, for the gene-level p in step 7.
# Why: step 6 tests every cis variant (median ~4,450 per gene) and step 7 multiplied the best p by that number
# (Bonferroni). Nearby variants are correlated (LD), so that over-penalises genes with dense LD; M_eff counts how
# many INDEPENDENT tests the gene's variants are worth (Li & Ji 2005: from the eigenvalues of their correlation).
# The same M_eff is used for the gene's real and shuffled scans (shuffling samples keeps the LD), so the
# calibration in step 7 stays valid. (caNRD_editQTL uses the same p x M_eff for edQTL sites.)
# Same variants as step 6: within 1 Mb of the TSS, MAF >= 0.05, the step-6 samples (all 670 donors).
# Li & Ji: M_eff = sum over eigenvalues lambda of [ (lambda >= 1) + (lambda - floor(lambda)) ]
# Eigenvalues from an SVD of the standardised genotypes (samples x variants), cheap since samples << variants.
# Inputs : output/06_scan_<CHR>_genes.csv (genes and their n_variants, as a check), output/genotypes/<CHR>, GENCODE
# Output : output/06b_meff_<CHR>.csv (gene_name, n_variants, m_eff)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))
library(parallel)

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHR      <- Sys.getenv("CHR", "chr22")
gtf_file <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf"
WINDOW   <- 1e6
MIN_MAF  <- 0.05
N_CORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

## ================= genes of the step-6 scan, with their TSS =================
scan <- fread(file.path(out_dir, paste0("06_scan_", CHR, "_genes.csv")))
scan <- unique(scan[set == "genotype", .(gene_name, n_variants)])
gtf <- fread(cmd = paste("grep -P '\\tgene\\t'", gtf_file), header = FALSE)
gtf$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9)
gtf$tss <- ifelse(gtf$V7 == "+", gtf$V4, gtf$V5)
genes <- gtf %>% filter(V1 == CHR, gene_name %in% scan$gene_name) %>% select(gene_name, tss)
genes <- genes[!duplicated(genes$gene_name), ]
nrow(genes)

## ================= genotypes: same filters as step 6 =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
geno_prefix <- file.path(out_dir, "genotypes", CHR)
pl <- read.plink(paste0(geno_prefix, ".bed"), paste0(geno_prefix, ".bim"), paste0(geno_prefix, ".fam"))
G_all <- as(pl$genotypes, "numeric")
G_all <- G_all[match(cov_df$SUBJID, rownames(G_all)), ]
af <- colMeans(G_all, na.rm = TRUE) / 2
keep <- pmin(af, 1 - af) >= MIN_MAF
G_all <- G_all[, keep]
pos <- pl$map$position[keep]
missing <- which(is.na(G_all), arr.ind = TRUE)
G_all[missing] <- colMeans(G_all, na.rm = TRUE)[missing[, 2]]
dim(G_all)
# (allele flipping does not change correlations, so it is not needed here)

## ================= M_eff per gene =================
li_ji <- function(lambda) sum((lambda >= 1) + (lambda - floor(lambda)))
meff_gene <- function(j) {
  cis <- which(abs(pos - genes$tss[j]) <= WINDOW)
  if (length(cis) == 0) return(NULL)
  G <- G_all[, cis, drop = FALSE]
  G <- G[, apply(G, 2, sd) > 0, drop = FALSE]
  Z <- scale(G)
  # eigenvalues of the variants' correlation matrix = squared singular values / (n - 1)
  lambda <- svd(Z, nu = 0, nv = 0)$d^2 / (nrow(Z) - 1)
  data.frame(gene_name = genes$gene_name[j], n_variants_6b = length(cis), m_eff = li_ji(lambda))
}
meff <- bind_rows(mclapply(seq_len(nrow(genes)), meff_gene, mc.cores = N_CORES, mc.preschedule = FALSE))
meff <- meff %>% left_join(scan, by = "gene_name")
# check: same variant count as step 6
cat("genes:", nrow(meff), "| variant count equals step 6:", mean(meff$n_variants_6b == meff$n_variants), "\n")
print(summary(meff$m_eff))
fwrite(meff, file.path(out_dir, paste0("06b_meff_", CHR, ".csv")))
message("06b_meff done: ", CHR)
