#author: haim krupkin+claude
#date: 11.13.2026
#description: step 10c. Find neutrophil eQTLs whose genotype lines visibly FAN OUT (as Kim-Hellmuth 2020 Fig. 1C), and
# plot the best 3 on the linear scale with the binned per-allele effect.
# Why neutrophils: their proportion varies the most between GTEx donors (most room for the lines to spread).
# Candidates (neutrophil eGenes, step 10b's candidate table):
#   - the gene is also made by other cell types (n_other_tested >= 1) and neutrophils make 30-80% of its RNA on
#     average (a gene made only by neutrophils gives a shift, not a fan)
#   - effect in neutrophils only (the other tested cell types |z| < 2) and replicated pre / post-mortem
#   - common variant: >= 30 donors in the smallest genotype group
# Score (TPM scale, x = neutrophil proportion, donors in 5 equal groups): the per-allele effect in the top fifth vs
# the bottom fifth (z of the difference); the bottom fifth should be near 0 (|effect| < 30% of the top fifth).
# Flag: also a Kim-Hellmuth 2020 neutrophil ieGene with the same lead signal (r2 >= 0.8; step 7c).
# Inputs : figures<SCAN_TAG>/variants/example_candidates_all_hits.csv (step 10b), 05b_mean_phi.csv,
#          07c_kimhellmuth_pairs<SCAN_TAG>.csv, and everything 10_plot_variant.R reads
# Outputs: figures<SCAN_TAG>/variants/10c_fan_candidates.csv (all scored), and the plots of the best 3
#          (<gene>_Neutrophils_<variant>_linear_binned.png)
# Run: Rscript 10c_pick_fan_examples.R
suppressMessages(library(dplyr))

## ================= load the plot functions (10_plot_variant.R) without running its own plots =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
Sys.setenv(PLOT_VARIANT_NO_MAIN = "1")
source(file.path(pipe_dir, "10_plot_variant.R"))
N_BEST    <- 3
MAX_TRIED <- 200          # score at most the 200 strongest candidates (each needs its genotypes read)

## ================= 1. candidates =================
cand <- read.csv(file.path(plot_dir, "example_candidates_all_hits.csv"))
cand <- cand %>%
  filter(celltype == "Neutrophils",
         n_other_tested >= 1,          # other cell types make the gene too
         max_other_z < 2,              # but the effect is in neutrophils only
         replicates)                   # same direction pre- and post-mortem

# the neutrophils' mean share of the gene's RNA
neut_share <- mean_phi %>% select(gene_name, neut_share = Neutrophils)
cand <- cand %>%
  left_join(neut_share, by = "gene_name") %>%
  filter(neut_share >= 0.3, neut_share <= 0.8)

# Kim-Hellmuth flag: their neutrophil ieGene, same lead signal as ours
kh <- read.csv(file.path(out_dir, paste0("07c_kimhellmuth_pairs", SCAN_TAG, ".csv")))
kh_genes <- kh %>%
  filter(their_iegene, our_egene, r2 >= 0.8) %>%
  pull(gene_name)
cand <- cand %>%
  mutate(kim_hellmuth_same_signal = gene_name %in% kh_genes) %>%
  arrange(desc(abs(z))) %>%
  slice_head(n = MAX_TRIED)
cat("candidates scored:", nrow(cand), "\n")

## ================= 2. score the fan of each candidate (TPM scale, x = neutrophil proportion) =================
scores <- list()
for (r in seq_len(nrow(cand))) {
  gene    <- cand$gene_name[r]
  variant <- cand$variant_id[r]
  res <- tryCatch({
    chr <- find_chr(gene, variant)
    d <- variant_data(gene, "Neutrophils", variant, chr, scale = "linear")
    d$x <- d$prop
    smallest_group <- min(table(d$genotype))
    b <- binned_effect(d, n_bins = 5)
    bottom <- b %>% filter(bin == 1)
    top    <- b %>% filter(bin == 5)
    data.frame(gene_name = gene,
               variant_id = variant,
               min_genotype_n = smallest_group,
               effect_bottom_fifth = bottom$effect,
               effect_top_fifth = top$effect,
               fan_z = (top$effect - bottom$effect) / sqrt(top$se^2 + bottom$se^2),
               bottom_vs_top = abs(bottom$effect) / abs(top$effect),
               steps_in_direction = sum(sign(diff(b$effect)) == sign(top$effect)))   # of 4 steps
  }, error = function(e) {
    message("skipped ", gene, ": ", conditionMessage(e))
    NULL
  })
  scores[[r]] <- res
}
scores <- bind_rows(scores)

## ================= 3. rank =================
ranked <- cand %>%
  select(gene_name, variant_id, z, max_other_z, n_other_tested, neut_share, kim_hellmuth_same_signal) %>%
  inner_join(scores, by = c("gene_name", "variant_id"))
ranked <- ranked %>%
  mutate(good_fan = min_genotype_n >= 30 &
                    sign(effect_top_fifth) == sign(z) &      # the fan points the way of the model's effect
                    bottom_vs_top < 0.3 &                    # ~0 where there are few neutrophils
                    steps_in_direction >= 3)                 # rises steadily (3 of 4 steps)
ranked <- ranked %>%
  arrange(desc(good_fan), desc(abs(fan_z)))
out_file <- file.path(plot_dir, "10c_fan_candidates.csv")
write.csv(ranked, out_file, row.names = FALSE)
cat("scored candidates:", out_file, "\n")
print(head(ranked, 10), digits = 2)

## ================= 4. plot the best 3: TPM scale + binned per-allele effect =================
best <- ranked %>%
  filter(good_fan) %>%
  slice_head(n = N_BEST)
# each: [NNLS expression by cell type] [TPM vs neutrophil proportion, per genotype + 95% CI] [binned per-allele effect]
for (r in seq_len(nrow(best))) {
  invisible(plot_gene(best$gene_name[r], best$variant_id[r], celltypes = "Neutrophils", profile = TRUE,
                      xaxis = "prop", scale = "linear", binned = TRUE))
}
