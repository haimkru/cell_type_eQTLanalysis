#author: haim krupkin+claude
#date: 11.10.2026
#description: step 8b. Power check (spike-in) results. For each spiked gene (step 8a plan, step 6 MODEL=spike scan):
#   detected (calibrated) : gene-level p of the spiked cell type's test is smaller than EVERY shuffled-genotype gene-level
#                           p of that test genome-wide (5 shuffles x all genes; the single-gene version of step 7's rule)
#   detected (nominal)    : gene-level p < 0.05 (Bonferroni within the gene only)
#   right cell type       : at that best variant, the spiked cell type has the largest |z| of the 6
#   leak                  : share of the OTHER 5 cell types with p < 0.05 at that variant
# beta_sd = 0 rows are the null control (should not be detected).
# Inputs : output/08_spike_plan_<CHR>_with_variants.csv, output/06_scan_<CHR>_{genes,top_variants}_spike.csv,
#          output/06_scan_chr*_genes.csv (main scans of all chromosomes: shuffled sets = the null)
# Output : output/08_spikein_power.csv (per cell type x effect size), output/08_spikein_genes.csv (per gene)
suppressMessages(library(dplyr))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHR      <- Sys.getenv("CHR", "chr22")
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")

# 11.10: NULL_TAG = which main scans give the null (e.g. "_phi_win100k_shuf25"); SPIKE_TAG = suffix of the spike scan
NULL_TAG  <- Sys.getenv("NULL_TAG", "_epc20")    # 11.13: main result = 20 hidden-noise PCs
SPIKE_TAG <- Sys.getenv("SPIKE_TAG", "_epc20")   # 11.13: main result = 20 hidden-noise PCs

## ================= the null: smallest shuffled gene-level p per test, genome-wide =================
main_files <- list.files(out_dir, pattern = paste0("^06_scan_chr[0-9]+_genes", NULL_TAG, "\\.csv$"), full.names = TRUE)
length(main_files)
null_min <- bind_rows(lapply(main_files, read.csv)) %>%
  filter(grepl("^shuffled", set), test %in% CELLTYPES) %>%
  group_by(test) %>%
  summarise(n_null = n(), p_null_min = min(p_gene), .groups = "drop")
print(as.data.frame(null_min))

## ================= spiked genes: result of the spiked cell type's test =================
plan  <- read.csv(file.path(out_dir, paste0("08_spike_plan_", CHR, "_with_variants", SPIKE_TAG, ".csv")))
genes <- read.csv(file.path(out_dir, paste0("06_scan_", CHR, "_genes_spike", SPIKE_TAG, ".csv")))
top   <- read.csv(file.path(out_dir, paste0("06_scan_", CHR, "_top_variants_spike", SPIKE_TAG, ".csv")))

res <- plan %>%
  left_join(top %>% filter(set == "genotype"), by = c("gene_name", "celltype" = "test"), suffix = c("_spiked", "")) %>%
  left_join(null_min, by = c("celltype" = "test"))
nrow(res)
sum(is.na(res$p_gene))          # genes with no result (should be 0)

# z of each cell type at the best variant of the spiked cell type's test
z <- as.matrix(res[, paste0("b_", CELLTYPES)]) / as.matrix(res[, paste0("se_", CELLTYPES)])
p <- as.matrix(res[, paste0("p_", CELLTYPES)])
target <- match(res$celltype, CELLTYPES)
res$detected_calibrated <- res$p_gene < res$p_null_min
res$detected_nominal    <- res$p_gene < 0.05
# (PHI: gated-out cell types are NA and are ignored)
res$right_celltype      <- apply(abs(z), 1, function(r) if (all(is.na(r))) NA else which.max(r)) == target
res$leak_other_p05      <- sapply(seq_len(nrow(res)), function(i) mean(p[i, -target[i]] < 0.05, na.rm = TRUE))
res$found_spiked_variant <- res$variant_id == res$variant_id_spiked
write.csv(res %>% select(gene_name, celltype, beta_sd, variant_id_spiked, variant_id, p_gene, p_null_min,
                         detected_calibrated, detected_nominal, right_celltype, leak_other_p05),
          file.path(out_dir, paste0("08_spikein_genes", SPIKE_TAG, ".csv")), row.names = FALSE)

## ================= power per cell type x effect size =================
power <- res %>%
  group_by(celltype, beta_sd) %>%
  # (new names: summarise would otherwise overwrite detected_nominal before the last two lines use it)
  summarise(genes = n(),
            share_detected_calibrated = mean(detected_calibrated),
            share_detected_nominal = mean(detected_nominal),
            n_detected_nominal = sum(detected_nominal),
            right_celltype_when_detected = if (any(detected_nominal)) mean(right_celltype[detected_nominal]) else NA,
            leak_other_p05_when_detected = if (any(detected_nominal)) mean(leak_other_p05[detected_nominal]) else NA,
            .groups = "drop") %>%
  arrange(factor(celltype, levels = CELLTYPES), beta_sd)
print(as.data.frame(power %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
write.csv(power, file.path(out_dir, paste0("08_spikein_power", SPIKE_TAG, ".csv")), row.names = FALSE)
