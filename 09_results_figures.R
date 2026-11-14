#author: haim krupkin+claude
#start date: 11.11.2026
#manual verification by hk: plots were manually verified starting on 11/12/2026
#code to generate data for plots wasnt yet manually verified
#description:we are testing if the eQTLs form this type of analysis worked and are robust.
#there is overal like 10 checks we need to do to ensure the hits are robust

#this methodology requires significant refrences to:
#validationL onek1k
#CSeQTL
#joint modeling (TCA)
#using proportions for cell type eQTLs-Westra 2015 and Zhernakova 2017,
#and the GTEx cell-type interaction eQTLs,
#Kim-Hellmuth 2020.
#permutation calibration - fastQTL


#   1. proportions      : the 6 cell-type proportions per sample (CIBERSORTx, batch-removed mixture)
#   2. batch            : variance of the proportions explained by each batch variable, raw vs batch-removed mixture;
#                         and the proportions pre- vs post-mortem. NOTE: way of death (+ centre, confounded with it) is
#                         removed by design in step 3, so its drop is not evidence; the corrected proportions are for the
#                         eQTL model only, not for studying pre- vs post-mortem composition (the gap may be real biology)
#   3. hits             : calibrated eGenes per cell type; the top hits and their effect (z) in every cell type
#   4. shuffle          : real genotypes vs shuffled genotypes (QQ plot, and eGenes called per held-out shuffle)
#   5. replication      : effect of each eGene in pre-mortem vs post-mortem draws, fitted separately
#   6. robustness       : effect with proportions from the batch-removed mixture vs from the raw mixture
#   7. OneK1K           : are our eGenes more often OneK1K eGenes (matching cell type) than the other genes?
#   8. spike-in         : power: known effects added to real expression (chr22), share found
#   9. OneK1K leads     : for eGenes also in OneK1K, is it the same signal (LD of the two lead variants) and direction (step 7b)
#  10. lab batch        : eGene effects with vs without the GTEx expression batch as covariate (chr21-22, LABBATCH=1 run)
#  11. Kim-Hellmuth     : our neutrophil eGenes vs the GTEx neutrophil interaction eQTLs (step 7c, if run)
#  12. ASE              : allelic ratio of heterozygotes vs the cell type's proportion (step 11, if run)
#  13. leakage          : r_b (true-effect correlation, Qi et al. 2018) between every pair of cell types, both
#                         directions; lineage structure; direct test on genes made almost only by one cell type
#  14. bulk overlap     : UpSet of the 6 cell types' eGenes and all GTEx bulk whole-blood eGenes
#  15. per-variant      : expression vs cell-type proportion by genotype for the genes in
#                         figures<SCAN_TAG>/variants/example_variants.csv (step 10b, or edit it by hand)
# eGene = calibrated as in step 7: each gene's p vs the pooled gene-level p of the 5 shuffled-genotype scans, then BH < 5%.
# Inputs : output/05_proportions_6classes{,_raw}.csv, 02_covariates.csv, GTEx sample attributes,
#          output/06_scan_<chr>_{genes,top_variants}{,_pre,_post,_rawprops}.csv, OneK1K, output/08_spikein_power.csv, output/07b_onek1k_lead_summary.csv
# Outputs: output/figures/09_<n>_<name>.png and .csv, and all figures in one PDF: output/figures/09_all_figures.pdf
#hk-first we load the packages
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ggpubr))      # stat_compare_means (figure 2b)

## ================= settings =================
#hk - these are general settings of variables
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
# 11.12: SCAN_TAG = the step-6 OUT_TAG to show (e.g. "_epc20"); figures then go to output/figures<SCAN_TAG>/
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")   # 11.13: main result = 20 hidden-noise PCs; "" = the model without
fig_dir  <- file.path(out_dir, paste0("figures", SCAN_TAG))
dir.create(fig_dir, showWarnings = FALSE)
ann_file <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/edQTLs_GTEx/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
#hk- to load the onek1k data
onek1k_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/onek1k_top_eqtl/unzipped/OneK1K_TensorQTL_top_eQTL_summary"

CHRS <- paste0("chr", 1:22)
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
FDR <- 0.05
#hk-this is to make the onek1k cell types consistent with our deconvulved cell types
#this needs to be done manually which kinda sad
ONEK1K_MATCH <- list(Bcells = c("B_IN", "B_MEM", "Plasma"), CD4 = c("CD4_NC", "CD4_ET", "CD4_SOX4"),
                     CD8 = c("CD8_NC", "CD8_ET", "CD8_S100B"), NK = c("NK", "NK_R"),
                     Monocytes = c("Mono_C", "Mono_NC"), Neutrophils = c("Mono_C"))
# one colour per cell type, used in every figure
#hk-we try to maintain consistency with the other plots we use
CT_COLS <- c(Neutrophils = "#1f77b4", Monocytes = "#ff7f0e", Bcells = "#2ca02c",
             CD4 = "#d62728", CD8 = "#9467bd", NK = "#8c564b")
#this is mostyl to get a consistent theme
theme_set(theme_bw(base_size = 12))

# save one figure (PNG) and the data behind it (CSV); print the full paths
#this is just a script to enable saving of pltos easily
#ideally we would also save .svg and .pdf but whatever
all_plots <- list()
save_fig <- function(p, data, name, width = 9, height = 6) {
  png_file <- file.path(fig_dir, paste0(name, ".png"))
  csv_file <- file.path(fig_dir, paste0(name, ".csv"))
  ggsave(png_file, p, width = width, height = height, dpi = 450)
  write.csv(data, csv_file, row.names = FALSE)
  cat("figure:", png_file, "\n  data:", csv_file, "\n")
  all_plots[[name]] <<- p
}

## ================= read the scans (all 22 chromosomes) =================
#hk- soem functions were ran chromsome wise, so we need to load each chromsome sepretly
read_scan <- function(kind, tag = "") {
  bind_rows(lapply(CHRS, function(ch) fread(file.path(out_dir, paste0("06_scan_", ch, "_", kind, tag, SCAN_TAG, ".csv")))))
}
#hk - this part can take a bit of time to load data
genes <- read_scan("genes")                     # one row per gene x test x set (real genotypes + 5 shuffles)
top   <- read_scan("top_variants")              # the best variant per gene x test x set, with all 6 effects
pre   <- read_scan("top_variants", "_pre")      # refit at the main best variants: pre-mortem draws only
post  <- read_scan("top_variants", "_post")     # ... post-mortem draws only
rawp  <- read_scan("top_variants", "_rawprops") # ... proportions from the raw (not batch-removed) mixture
genes <- genes %>% filter(test %in% CELLTYPES)
top   <- top %>% filter(set == "genotype", test %in% CELLTYPES)
cat("genes per set:\n")
print(genes %>% distinct(gene_name, set) %>% count(set))

## ================= calibrated eGenes (as step 7) =================
# empirical p of a gene = (1 + number of shuffled gene-level p <= its p) / (1 + number of shuffled p); then BH
#hk-to get the empirical fdr FOR EACH test by testing against ashuffled gene level variant, so its a 
#pertumation based FDR.
empirical_fdr <- function(p_obs, p_null) {
  p_null <- sort(p_null)
  p_emp <- (1 + findInterval(p_obs, p_null)) / (1 + length(p_null))
  p.adjust(p_emp, "BH")
}
genes
#HK-we ran the analysis 6 times overall, so 1 with the true genotype and 5 with shuffled genotype, as our
#null model.

null_p <- genes %>% filter(set != "genotype") %>% select(test, p_gene)    # null = all shuffles of this test
#variant level
real <- genes %>%
  filter(set == "genotype") %>%
  group_by(test) %>%
  mutate(fdr_emp = empirical_fdr(p_gene, null_p$p_gene[null_p$test == first(test)])) %>%
  ungroup()
#this dta frame saves the variant levle hits
real
#HK- this is us getting the number of hits after using a permutation based apporach to get p.values as well as to 
#adjust for multiple testing by benjamini hochberg.
#this part is a part that we developed here-  so we should sit and think about it a bit.
#this methodology requires signficant refrences.
#ok so the p.value correction is done by pertumation per gene - as many variants are in LD
#and the BH adjust for multiple testing (as there are many genes)

#eGene level hits, so we filter by FDR after doing the whole empircal p.value from permutation
#
hits <- real %>%
  filter(fdr_emp < FDR) %>%
  select(chr, gene_name, celltype = test, p_gene, fdr_emp)
hits$celltype <- factor(hits$celltype, levels = CELLTYPES)
print(table(hits$celltype))

## ================= 1. proportions =================
props_main <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
props_raw  <- read.csv(file.path(out_dir, "05_proportions_6classes_raw.csv"))
props_long <- bind_rows(
  props_main %>% pivot_longer(all_of(CELLTYPES), names_to = "celltype", values_to = "proportion") %>%
    mutate(mixture = "batch removed (used)"),
  props_raw %>% pivot_longer(all_of(CELLTYPES), names_to = "celltype", values_to = "proportion") %>%
    mutate(mixture = "raw (check)"))
props_long$celltype <- factor(props_long$celltype, levels = CELLTYPES)

# 1a: distribution per cell type
p1a <- ggplot(props_long %>% filter(mixture == "batch removed (used)"), aes(celltype, proportion, fill = celltype)) +
  geom_violin(scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  scale_fill_manual(values = CT_COLS, guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "proportion per sample", title = "Cell-type proportions (CIBERSORTx, 670 GTEx blood samples)")
#these are the cell type proportions after batch effect correction, 
#i think they are close to the real range
#and thus i trust them
p1a
# 1b: every sample as one stacked bar, ordered by neutrophil share
sample_order <- props_main$SAMPID[order(props_main$Neutrophils)]
#hk - this plot isnt very important ngl
p1b <- ggplot(props_long %>% filter(mixture == "batch removed (used)") %>%
                mutate(SAMPID = factor(SAMPID, levels = sample_order)),
              aes(SAMPID, proportion, fill = celltype)) +
  geom_col(width = 1) +
  scale_fill_manual(values = CT_COLS) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "samples (ordered by neutrophil share)", y = NULL, fill = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p1b
prop_summary <- props_long %>%
  group_by(mixture, celltype) %>%
  summarise(mean = mean(proportion), median = median(proportion),
            q05 = quantile(proportion, 0.05), q95 = quantile(proportion, 0.95), .groups = "drop")
print(as.data.frame(prop_summary), digits = 2)
save_fig(p1a / p1b, prop_summary, "09_1_proportions", width = 10, height = 8)
save_fig(p1a / p1b, prop_summary, "09_1_proportions", width = 10, height = 8)

## ================= 2. batch effects in the proportions =================
# batch variables: way of death (pre/post-mortem draw, ventilator) and lab batches (centre, nucleic-acid isolation
# batch SMNABTCH, expression batch SMGEBTCH)
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv")) %>% select(SAMPID, post_mortem, ventilator)
ann <- fread(ann_file, select = c("SAMPID", "SMCENTER", "SMNABTCH", "SMGEBTCH"), quote = "")
batch_df <- cov_df %>%
  left_join(ann, by = "SAMPID") %>%
  mutate(post_mortem = ifelse(post_mortem == 1, "post-mortem", "pre-mortem"),
         ventilator = ifelse(round(ventilator) == 1, "ventilator", "other death"))   # missing ventilator -> filled mean
BATCH_VARS <- c("post_mortem", "ventilator", "SMCENTER", "SMNABTCH", "SMGEBTCH")

# variance of each proportion explained by each batch variable:
#adjusted R2 (adjusted, because the lab batches have
# many levels and plain R2 grows with the number of levels even without an effect)
#this is us checking if we were able to remove the batch effects.
#which seems like we were largely able to make huge white bars dissapear
batch_r2 <- props_long %>%
  left_join(batch_df, by = "SAMPID") %>%
  pivot_longer(all_of(BATCH_VARS), names_to = "batch_variable", values_to = "batch") %>%
  group_by(mixture, celltype, batch_variable) %>%
  summarise(adj_r2 = summary(lm(proportion ~ factor(batch)))$adj.r.squared,
            n_levels = n_distinct(batch), .groups = "drop") %>%
  mutate(batch_variable = factor(batch_variable, levels = BATCH_VARS))
p2a <- ggplot(batch_r2, aes(batch_variable,
                            pmax(adj_r2, 0),
                            fill = mixture)) +
  geom_col(position = "dodge") +
  facet_wrap(~ celltype, nrow = 1) +
  scale_fill_manual(values = c("batch removed (used)" = "grey30", "raw (check)" = "grey75")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample Metadata", y = "variance explained (adj. R2)", fill = "Analysis_Method") +
  theme(plot.subtitle = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
#hk- this plot is pretty important andi do want to save it :)
p2a
p2a
# the largest batch effect (way of death): proportions pre- vs post-mortem, raw vs batch-removed
plot_df <- props_long %>%
  left_join(batch_df %>% select(SAMPID, post_mortem), by = "SAMPID") %>%
  filter(!is.na(post_mortem))

# t-test per cell type x mixture, done by hand
stars_df <- plot_df %>%
  group_by(mixture, celltype) %>%
  summarise(p = t.test(proportion ~ post_mortem)$p.value, .groups = "drop") %>%
  mutate(label = as.character(cut(p, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                                  labels = c("****", "***", "**", "*", "ns"))))

p2b <- ggplot(plot_df, aes(celltype, proportion, fill = post_mortem)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_text(data = stars_df, aes(celltype, 1.06, label = label),
            inherit.aes = FALSE, size = 4) +
  facet_wrap(~ mixture) +
  scale_fill_manual(values = c("pre-mortem" = "#4daf4a", "post-mortem" = "#984ea3")) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     labels = scales::percent) +
  labs(x = NULL, y = "Cell Type Proportion in Sample", fill = "Blood_Draw_Timing") +
  coord_flip() +                                   # no clip = "off"
  theme(legend.position = "bottom",
        plot.margin = margin(10, 20, 10, 10))
p2b

save_fig(p2a, stars_df, "PRE_VERSUS_POST_varaince_explained_by_each_variable_after_batch_correction",
         width = 8, height = 4)
save_fig(p2b, stars_df, "pre_versus_post_batch_effect_correction_cell_type_proportions",
         width = 5, height = 4)

#save_fig(p2a / p2b, batch_r2, "09_2_batch_effects", width = 14, height = 9)

## ================= 3. the hits =================
# 3a: number of calibrated eGenes per cell type
n_hits <- hits %>% count(celltype, name = "egenes", .drop = FALSE) %>%
  left_join(real %>% count(celltype = factor(test, levels = CELLTYPES), name = "genes_tested"), by = "celltype")
p3a <- ggplot(n_hits, aes(celltype, egenes, fill = celltype)) +
  geom_col() +
  geom_text(aes(label = egenes), vjust = -0.3) +
  scale_fill_manual(values = CT_COLS, guide = "none") +
  scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "N eGenes (FDR < 5%)")
p3a
save_fig(p3a, stars_df, "per_cel_type_number_of_eGenes",
         width = 6, height = 4)
# 3b: the top 8 eGenes of each cell type: z (effect / SE)
#of every cell type at the gene's best variant.
#so we sellect e genes by theire z score so their effect compared to their standard error
# blank = cell type not tested in that gene (it carries < 10% of the gene's RNA)
top_hits <- hits %>% group_by(celltype) %>% arrange(p_gene) %>% slice_head(n = 5) %>% ungroup()
z_long <- top_hits %>%
  left_join(top %>% mutate(celltype = factor(test, levels = CELLTYPES)) %>% select(-p_gene),
            by = c("gene_name", "celltype")) %>%
  select(called_in = celltype, gene_name, variant_id, starts_with("b_"), starts_with("se_")) %>%
  pivot_longer(-c(called_in, gene_name, variant_id), names_to = c(".value", "effect_in"), names_sep = "_") %>%
  mutate(z = b / se,
         effect_in = factor(effect_in, levels = CELLTYPES),
         gene_label = factor(paste0(gene_name, " (", called_in, ")"),
                             levels = rev(unique(paste0(gene_name, " (", called_in, ")")))))
p3b <- ggplot(z_long, aes(effect_in, gene_label, fill = pmax(pmin(z, 10), -10))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(is.na(z), "", round(z, 1))), size = 2.6) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", na.value = "grey90",
                       limits = c(-10, 10), name = "z") +
  labs(x = "Cell Type eGene", y = NULL       ) +
  theme(panel.grid = element_blank())
p3b
save_fig(p3b, stars_df, "Cell_Type_eGene",
         width = 8, height = 4)
print(as.data.frame(n_hits))

hits <- hits %>%
  mutate(receptor_gene = case_when(grepl("^IG[HKL][VDJC]|^IGLL|^IGH(M|D|E|G[1-4]|A[12])$", gene_name) ~ "antibody (B)",
                                   grepl("^TR[ABDG][VDJC]", gene_name) ~ "T-cell receptor (T)",
                                   TRUE ~ ""),
         lineage_conflict = (receptor_gene == "antibody (B)" & celltype != "Bcells") |
           (receptor_gene == "T-cell receptor (T)" & !celltype %in% c("CD4", "CD8")))
cat("receptor-gene hits per cell type, and lineage conflicts:\n")
print(hits %>% filter(receptor_gene != "") %>% count(celltype, receptor_gene, lineage_conflict) %>% as.data.frame())
write.csv(hits %>% arrange(celltype, p_gene), file.path(fig_dir, "09_3_all_hits.csv"), row.names = FALSE)
cat("all hits:", file.path(fig_dir, "09_3_all_hits.csv"), "\n")
save_fig(p3a | p3b, z_long, "09_3_hits", width = 15, height = 11)

## ================= 4. real vs shuffled genotypes =================
# 4a: QQ plot of the gene-level p: real genotypes vs each shuffle (shuffles should follow the diagonal)
qq <- genes %>%
  mutate(celltype = factor(test, levels = CELLTYPES),
         data = ifelse(set == "genotype", "real genotypes", "shuffled genotypes")) %>%
  group_by(celltype, set, data) %>%
  arrange(p_gene) %>%
  mutate(expected = -log10(ppoints(n())), observed = -log10(p_gene)) %>%
  ungroup()
qq_plot_data <- qq %>% filter(observed > 0.5 | row_number() %% 20 == 0)   # thin the dense bulk near 0 for plotting
p4a <- ggplot(qq_plot_data, aes(expected, observed, colour = data, group = set)) +
  geom_abline(linetype = "dashed") +
  geom_point(size = 0.5) +
  facet_wrap(~ celltype, nrow = 2) +
  scale_colour_manual(values = c("real genotypes" = "black", "shuffled genotypes" = "grey60")) +
  labs(x = "Expected -log10 p", y = "Observed -log10 p (gene level)", colour = NULL,) +
  theme(legend.position = "top")
p4a
###
#hk - we are making a slightly better more refined qqplot
# split real vs shuffled; shuffles are drawn first (grey) so the real points sit on top
real_qq <- qq_plot_data %>% filter(data == "real genotypes")
null_qq <- qq_plot_data %>% filter(data == "shuffled genotypes")

# the rightmost point of each cell type (largest expected -log10 p = its top gene)
label_df <- real_qq %>%
  group_by(celltype) %>%
  slice_max(expected, n = 1, with_ties = FALSE) %>%
  ungroup()

real_qq <- qq_plot_data %>% filter(data == "real genotypes")
null_qq <- qq_plot_data %>% filter(data == "shuffled genotypes")

# labels: the rightmost point of each cell type, plus the rightmost shuffled point
label_df <- bind_rows(
  real_qq %>%
    group_by(celltype) %>%
    slice_max(expected, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = as.character(celltype)),
  null_qq %>%
    slice_max(expected, n = 1, with_ties = FALSE) %>%
    mutate(label = "shuffled")
)

p4a <- ggplot() +
  geom_abline(linetype = "dashed") +
  geom_point(data = null_qq, aes(expected, observed), colour = "grey75", size = 0.5) +
  geom_point(data = real_qq, aes(expected, observed, colour = celltype), size = 0.6) +
  ggrepel::geom_text_repel(data = label_df,
                           aes(expected, observed, label = label, colour = label),
                           nudge_x = 0.3, hjust = 1, direction = "y",
                           segment.colour = NA, size = 3.5, fontface = "bold") +
  scale_colour_manual(values = c(CT_COLS, shuffled = "grey50"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = "Expected -log10 p", y = "Observed -log10 p (gene level)")
p4a
save_fig(p4a, stars_df, "qqplot_per_cell_type",
         width = 4, height = 4)
# 4b: eGenes called with real genotypes vs with each shuffle treated as real (scored against the other 4 shuffles)
shuffle_sets <- sort(unique(genes$set[genes$set != "genotype"]))
calls <- bind_rows(lapply(CELLTYPES, function(ct) {
  d <- genes %>% filter(test == ct)
  held_out <- sapply(shuffle_sets, function(s) {
    sum(empirical_fdr(d$p_gene[d$set == s], d$p_gene[d$set %in% setdiff(shuffle_sets, s)]) < FDR)
  })
  data.frame(celltype = ct, data = c("real", shuffle_sets),
             egenes = c(sum(hits$celltype == ct), held_out))
}))
calls$celltype <- factor(calls$celltype, levels = CELLTYPES)
p4b <- ggplot(calls, aes(data, egenes, fill = data == "real")) +
  geom_col() +
  geom_text(aes(label = egenes), vjust = -0.3, size = 3) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "calibrated eGenes", title = "eGenes: real genotypes vs each held-out shuffle") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4b
print(calls %>% pivot_wider(names_from = data, values_from = egenes) %>% as.data.frame())
save_fig(p4a / p4b + plot_layout(heights = c(2, 1)), calls, "09_4_real_vs_shuffled", width = 13, height = 11)
p4b
## ================= 5. replication: pre- vs post-mortem draws =================
#here is just further replication that indeed, the effect size 
# effect of the called cell type at the eGene's best variant, fitted separately in each group of draws
effect_at <- function(d, label) {
  d %>% filter(set == "genotype") %>%
    distinct(gene_name, variant_id, .keep_all = TRUE) %>%
    select(gene_name, variant_id, all_of(paste0("b_", CELLTYPES))) %>%
    pivot_longer(-c(gene_name, variant_id), names_to = "celltype", values_to = label) %>%
    mutate(celltype = sub("^b_", "", celltype))
}
hit_effects <- hits %>%
  mutate(celltype = as.character(celltype)) %>%
  left_join(top %>% select(gene_name, celltype = test, variant_id), by = c("gene_name", "celltype")) %>%
  left_join(effect_at(top, "effect_main"), by = c("gene_name", "variant_id", "celltype")) %>%
  left_join(effect_at(pre, "effect_pre"), by = c("gene_name", "variant_id", "celltype")) %>%
  left_join(effect_at(post, "effect_post"), by = c("gene_name", "variant_id", "celltype")) %>%
  left_join(effect_at(rawp, "effect_raw"), by = c("gene_name", "variant_id", "celltype")) %>%
  mutate(celltype = factor(celltype, levels = CELLTYPES))
# one label per panel: genes with both fits, share with the same sign, Spearman rho
panel_stats <- function(d, x, y) {
  d %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    group_by(celltype) %>%
    summarise(n = n(), same_sign = mean(sign(.data[[x]]) == sign(.data[[y]])),
              rho = cor(.data[[x]], .data[[y]], method = "spearman"), .groups = "drop") %>%
    mutate(label = sprintf("n = %d\nsame sign %.0f%%\nSpearman Rho = %.2f", n, 100 * same_sign, rho))
}
rep_stats <- panel_stats(hit_effects, "effect_pre", "effect_post")
# for display only: a few fits in the smaller subsets are unstable (huge effects), they are drawn at the edge (+-4);
# the statistics above use the real values (rho is rank-based, so those few points do not drive it)
CLIP <- 4
clip <- function(x) pmax(pmin(x, CLIP), -CLIP)
p5 <- ggplot(hit_effects %>% filter(!is.na(effect_pre), !is.na(effect_post)),
             aes(clip(effect_pre), clip(effect_post), colour = celltype)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(linetype = "dashed") +
  geom_point(size = 1, alpha = 0.6) +
  geom_text(data = rep_stats, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.1,
            colour = "black", size = 2) +
  facet_wrap(~ celltype, nrow = 2) +
  scale_colour_manual(values = CT_COLS, guide = "none") +
  labs(x = "effect, pre-mortem draws only", y = "effect, post-mortem draws only")
p5
print(as.data.frame(rep_stats %>% select(-label)), digits = 2)
save_fig(p5, hit_effects, "effects_are_pre_or_post_only", width = 6, height = 4)

## ================= 6. robustness: raw-mixture proportions =================
#this isint that important as the previous plot already does the same job, so i would hide it.

rob_stats <- panel_stats(hit_effects, "effect_main", "effect_raw")
p6 <- ggplot(hit_effects %>% filter(!is.na(effect_raw)), aes(clip(effect_main), clip(effect_raw), colour = celltype)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(linetype = "dashed") +
  geom_point(size = 1, alpha = 0.6) +
  geom_text(data = rob_stats, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.1,
            colour = "black", size = 3) +
  facet_wrap(~ celltype, nrow = 2) +
  scale_colour_manual(values = CT_COLS, guide = "none") +
  labs(x = "effect, proportions from the batch-removed mixture (main)",
       y = "effect, proportions from the raw mixture",
       title = "Robustness: eGene effects with the two sets of proportions",
       subtitle = "effects beyond +-4 drawn at the edge")
print(as.data.frame(rob_stats %>% select(-label)), digits = 2)
p6
save_fig(p6, rob_stats %>% select(-label), "09_6_robustness_raw_proportions", width = 11, height = 7)

## ================= 7. OneK1K gene overlap =================
#this is super imporatn as its a validaiton figure
# OneK1K eGenes (q < 0.05) in the matching cell types; neutrophils (not in OneK1K, a PBMC study) vs classical monocytes
# OneK1K eGenes (q < 0.05) in the matching cell types;
#neutrophils (not in OneK1K, a PBMC study) vs classical monocytes
onek1k <- bind_rows(lapply(unique(unlist(ONEK1K_MATCH)), function(ty) {
  bind_rows(lapply(CHRS, function(ch) {
    fread(file.path(onek1k_dir, paste0("OneK1K_", ty, ".sig_cis_qtl_pairs.", ch, ".csv")),
          select = c("phenotype_id", "qval")) %>%
      filter(qval < FDR) %>%
      distinct(phenotype_id) %>%
      mutate(onek1k_type = ty)
  }))
}))

ok_tab <- bind_rows(lapply(CELLTYPES, function(ct) {
  ok_genes <- unique(onek1k$phenotype_id[onek1k$onek1k_type %in% ONEK1K_MATCH[[ct]]])
  d <- real %>%
    filter(test == ct) %>%
    mutate(our_egene = fdr_emp < FDR, in_onek1k = gene_name %in% ok_genes)
  # Fisher test once per cell type; factor levels force a full 2x2 table
  ft <- fisher.test(table(factor(d$our_egene, levels = c(FALSE, TRUE)),
                          factor(d$in_onek1k, levels = c(FALSE, TRUE))))
  d %>%
    group_by(our_egene) %>%
    summarise(genes = n(), share_in_onek1k = mean(in_onek1k), .groups = "drop") %>%
    mutate(celltype = ct,
           fisher_p = ft$p.value,
           odds_ratio = unname(ft$estimate))   # odds of being a OneK1K eGene: our eGenes vs other tested genes
}))
ok_tab <- ok_tab %>%
  mutate(group = ifelse(our_egene, "our eGenes", "other tested genes"),
         celltype = factor(celltype, levels = CELLTYPES))

# label: OR on top, p below; placed above the taller of the two bars
ok_p <- ok_tab %>%
  group_by(celltype) %>%
  summarise(y = max(share_in_onek1k),
            odds_ratio = first(odds_ratio),
            fisher_p = first(fisher_p), .groups = "drop") %>%
  mutate(p_label = ifelse(fisher_p < 0.005, "p < 0.005", paste0("p = ", signif(fisher_p, 2))),
         label = paste0("OR = ", sprintf("%.1f", odds_ratio), "\n", p_label))

p7 <- ggplot(ok_tab, aes(celltype, share_in_onek1k, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(data = ok_p, aes(x = celltype, y = y, label = label), vjust = -0.3, size = 3.2,
            lineheight = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c("our eGenes" = "black", "other tested genes" = "grey70")) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +  # room for 2-line labels
  labs(x = NULL, y = "eGene share overlap with OneK1K", fill = "Gene eGene or Not")+
  theme(legend.position = "bottom")
print(as.data.frame(ok_tab %>% select(celltype, group, genes, share_in_onek1k, odds_ratio, fisher_p)), digits = 2)
p7
save_fig(p7, ok_tab, "validation_with_onek1k1", width = 6, height = 5)
## ================= 8. spike-in power (chr22) =================
#this iis mainly to check leakage, which is the key thing for us!!

# 11.13: the spike-in simulation (step 8, chr22), refitted at the PLANTED variant (no variant search, so no selection):
# genes with no signal got a known effect (0.5-4 SD per allele) in ONE cell type that makes >= 20% of their RNA,
# through its RNA share (phi x genotype x effect), exactly as step 6 MODEL=spike. The joint model is refitted for
# each gene (10_plot_variant.R joint_refit: same design, weights and HC3 SEs as step 6, checked to give step 6's p).
#   8a: one planted gene, all its tested cell types: is the effect found only where it was planted?
#   8b: all planted genes: estimated vs true effect, for the planted cell type and for every other tested cell type
# NOTE: planted with the same proportions / RNA shares the model uses -> tests the model's own cross-talk, not errors
# of the deconvolution.
spike_plan_file <- file.path(out_dir, paste0("08_spike_plan_chr22_with_variants", SCAN_TAG, ".csv"))
if (file.exists(spike_plan_file)) {

  # load the refit function from step 10 (without running its own plots)
  Sys.setenv(PLOT_VARIANT_NO_MAIN = "1")
  pv8 <- new.env()
  sys.source(file.path(pipe_dir, "10_plot_variant.R"), envir = pv8)
  Sys.unsetenv("PLOT_VARIANT_NO_MAIN")

  # the plan: gene, planted cell type, effect size (SD per allele), planted variant
  spike_plan <- read.csv(spike_plan_file)
  spike_plan <- spike_plan %>% filter(beta_sd > 0)
  pv8$prefetch_dosages(spike_plan$variant_id, "chr22")      # read all planted variants at once

  # refit each planted gene at its planted variant; keep all 6 estimates
  spike_list <- list()
  for (r in seq_len(nrow(spike_plan))) {
    gene    <- spike_plan$gene_name[r]
    planted <- spike_plan$celltype[r]
    fit <- pv8$joint_refit(gene, spike_plan$variant_id[r], "chr22",
                           plant_celltype = planted, plant_sd = spike_plan$beta_sd[r])
    for (ct in CELLTYPES[fit$active]) {                        # only the cell types tested in this gene
      term <- paste0("g_", ct)
      spike_list[[length(spike_list) + 1]] <- data.frame(
        gene_name = gene, planted_in = planted, beta_sd = spike_plan$beta_sd[r],
        celltype = ct, role = ifelse(ct == planted, "Injected cell type", "Other cell types"),
        true_effect = ifelse(ct == planted, fit$true_effect, 0),
        planted_true_effect = fit$true_effect,
        estimate = fit$b[[term]], se = fit$se[[term]], p = fit$p[[term]])
    }
  }
  spike_fit <- bind_rows(spike_list)
  spike_fit <- spike_fit %>%
    mutate(role = factor(role, levels = c("Injected cell type", "Other cell types")))

  # 8b: all planted genes - estimated vs true effect; slope through 0 for each group
  # (planted cell type: on the diagonal = recovered; other cell types: flat at 0 = no leakage)
  slopes <- spike_fit %>%
    group_by(role) %>%
    summarise(slope = coef(lm(estimate ~ 0 + planted_true_effect))[1], .groups = "drop")
  slope_planted <- slopes$slope[slopes$role == "Injected cell type"]
  slope_others  <- slopes$slope[slopes$role == "Other cell types"]
  slope_text <- paste0("recovered: slope ", sprintf("%.2f", slope_planted), "\n",
                       "leaked: slope ", sprintf("%.2f", slope_others))
  cat("8b:", gsub("\n", "; ", slope_text), "\n")

  # same look as the other figures: black = ours / the signal, grey70 = the comparison; legend at the bottom
  p8b <- ggplot(spike_fit, aes(planted_true_effect, estimate, colour = role)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +   # perfect recovery
    geom_hline(yintercept = 0, colour = "grey50") +                                  # no leakage
    geom_point(size = 1, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = TRUE, linewidth = 1) +
    annotate("text", x = -Inf, y = Inf, label = slope_text, hjust = -0.05,
             vjust = 1.2, size = 5,
             lineheight = 0.9) +
    scale_colour_manual(values = c("Injected cell type" = "black", "Other cell types" = "grey70")) +
    labs(x = "True Injected Beta", y = "Estimated Beta", colour = NULL) +
    theme_bw(base_size = 12) +                     # set here too, so it looks the same when only section 8 is run
    theme(legend.position = "bottom")
  p8b
  save_fig(p8b, spike_fit, "spikein_leakage_recovered_vs_leaked", width = 4, height = 4)
}


## ================= 11. GTEx neutrophil interaction eQTLs (Kim-Hellmuth 2020; step 7c) =================
#only important part is 11a and should be in same foramt as p9
## ================= 11. GTEx neutrophil interaction eQTLs (Kim-Hellmuth 2020; step 7c) =================
kh_file <- file.path(out_dir, paste0("07c_kimhellmuth_summary", SCAN_TAG, ".csv"))
if (file.exists(kh_file)) {
  kh <- read.csv(kh_file)
  
  # 11a: share that are their neutrophil ieGenes: our neutrophil eGenes vs our other tested genes (same format as 7)
  kh_overlap <- data.frame(celltype = factor("Neutrophils", levels = CELLTYPES),
                           group = factor(c("our eGenes", "other tested genes"),
                                          levels = c("our eGenes", "other tested genes")),
                           share = c(kh$share_of_ours_in_theirs, kh$share_of_other_genes_in_theirs))
  kh_p <- data.frame(celltype = factor("Neutrophils", levels = CELLTYPES),
                     y = max(kh_overlap$share),
                     label = paste0("OR = ", sprintf("%.1f", kh$odds_ratio), "\n",
                                    ifelse(kh$fisher_p < 0.005, "p < 0.005",
                                           paste0("p = ", signif(kh$fisher_p, 2)))))
  
  p11a <- ggplot(kh_overlap, aes(celltype, share, fill = group)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(data = kh_p, aes(x = celltype, y = y, label = label), vjust = -0.3, size = 3.2,
              lineheight = 0.9, inherit.aes = FALSE) +
    scale_fill_manual(values = c("our eGenes" = "black", "other tested genes" = "grey70")) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "eGene share overlap with Kim-Hellmuth et al", fill = "Gene eGene or Not") +
    theme(legend.position = "bottom",
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 24))
  
  p11a
  save_fig(p11a, kh, "overlap_with_kim_hellmuth", width = 5, height = 6)
  
  # 11b: unchanged (kh_ld / p11b code as before), still saved with 11a:
  # save_fig(p11a | p11b, kh, "09_11_kimhellmuth_neutrophils", width = 12, height = 6)
}

## ================= 12. allele-specific expression (step 11) =================
#this is further valdiaiton but more like internal validaiton by using ASE which is pretty cool ngl

# drawn only if step 11 has run for this scan (step 11 needs this script's hit list, so run 9 -> 11 -> 9 again)
ase_file <- file.path(out_dir, paste0("11_ase_validation_summary", SCAN_TAG, ".csv"))
if (file.exists(ase_file)) {
  # read the step-11 summary (one row per cell type)
  ase_raw <- read.csv(ase_file)
  ase_raw <- ase_raw %>%
    mutate(celltype = factor(celltype, levels = CELLTYPES))
  
  # DENOMINATOR = ALL eGenes of the cell type, as in the other validation plots (e.g. OneK1K).
  # eGenes that ASE cannot test (too few heterozygotes or reads) count as "not validated".
  #   egenes          : all our eGenes in this cell type
  #   egenes_testable : the eGenes ASE could test
  #   n_validated     : testable eGenes validated by ASE (real proportions)
  #   share_validated_shuffled : share of the testable eGenes validated when the proportions are shuffled (chance)
  
  # number of validated eGenes: real, and chance (the chance share turned back into a count)
  ase_raw <- ase_raw %>%
    mutate(k_real = n_validated,
           k_shuf = round(share_validated_shuffled * egenes_testable))
  
  # the shares out of ALL eGenes (these are the bar heights)
  ase_raw <- ase_raw %>%
    mutate(share_real_of_egenes = k_real / egenes,
           share_shuf_of_egenes = k_shuf / egenes)
  
  # Fisher test per cell type: validated vs not validated, our eGenes vs chance (same denominator: all eGenes)
  fisher_or_p <- function(k1, k2, n) {
    if (is.na(n) || n == 0) return(c(NA, NA))
    ft <- fisher.test(matrix(c(k1, n - k1, k2, n - k2), nrow = 2, byrow = TRUE))
    c(unname(ft$estimate), ft$p.value)
  }
  res <- t(mapply(fisher_or_p, ase_raw$k_real, ase_raw$k_shuf, ase_raw$egenes))
  ase_raw$odds_ratio <- res[, 1]
  ase_raw$fisher_p   <- res[, 2]
  
  # long format for the bars: one row per cell type x group
  ase <- ase_raw %>%
    select(celltype, egenes,
           `our eGenes` = share_real_of_egenes,
           `chance` = share_shuf_of_egenes)
  ase <- ase %>%
    pivot_longer(c(`our eGenes`, `chance`), names_to = "group", values_to = "share")
  ase <- ase %>%
    mutate(group = factor(group, levels = c("our eGenes", "chance")))
  ase
  
  # label above the taller bar of each pair: "validated / all eGenes", then OR and p
  ase_p <- ase_raw %>%
    mutate(y = pmax(share_real_of_egenes, share_shuf_of_egenes),
           count_label = paste0(k_real, " / ", egenes),
           p_label = ifelse(fisher_p < 0.005, "p < 0.005", paste0("p = ", signif(fisher_p, 2))))
  ase_p <- ase_p %>%
    mutate(label = ifelse(is.na(fisher_p),
                          count_label,
                          paste0("OR = ", sprintf("%.1f", odds_ratio), "\n", p_label)))
  
  p12 <- ggplot(ase, aes(celltype, share, fill = group)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(data = ase_p, aes(x = celltype, y = y, label = label), vjust = -0.3, size = 3,
              lineheight = 0.9, inherit.aes = FALSE) +
    scale_fill_manual(values = c("our eGenes" = "black", "chance" = "grey70")) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.3))) +  # room for 3-line labels
    labs(x = NULL, y = "share of all eGenes validated via ASE", fill = "eGene or chance") +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))
  p12
}

## ================= 13. leakage: do effects cross between cell types? (step 7d) =================
#here we are checking for leakage

# the leaky CIBERSORTx/TCA profiles gave perviously 0.96-0.99 for every pair.
#Drawn only if step 7d has run for this scan.
# The 6 effects are estimated together from correlated proportions, so the estimation alone gives a negative
# correlation (shown by the shuffled scans); real biology: related cell types (T/NK; neutrophils-monocytes) share more.
#
# 11.13: the check is now r_b (Qi et al. 2018, Nat Commun 9:2282) instead of a plain correlation / |z| > 2 threshold.
# r_b = the correlation of the TRUE effects of two cell types, corrected for
#   (1) the estimation error of each effect (its SE): plain correlations are pulled toward 0 by noise
#   (2) correlated errors: all 6 effects come from one model on the same donors, so their errors are correlated
#       (re, the "sample overlap" term of r_b). re is taken from the shuffled-genotype scans (no true effects there).
#   r_b = [cov(b1, b2) - re * mean(s1 * s2)] / sqrt([var(b1) - mean(s1^2)] * [var(b2) - mean(s2^2)])
# Expected without leakage: related cell types (CD4 / CD8 / NK) high, myeloid vs lymphoid lower. Most cis-eQTLs are
# shared, so distant pairs are not expected near 0; leakage = a FLAT matrix (all pairs alike, as the old profiles).
# Limits: (a) genes are chosen because they are significant in the discovery cell type (winner's curse): its effects
#   are inflated, which pulls r_b toward 0. Both directions are shown (discover in A test in B, and the reverse);
#   if they agree, the bias is small. (b) re is one value per pair, pooled over genes. (c) Pairs with few genes are noisy.

## ================= 14. discovery overlap with bulk: UpSet =================
# all GTEx bulk whole-blood eGenes (qval <= 0.05) among the tested genes, and our eGenes of every cell type;
# every bar labelled with its count (the ~9,300 bulk-only genes would otherwise hide the small groups)
bulk_file <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz"

# the bulk eGenes, only those we also tested
bulk_table <- fread(bulk_file, select = c("gene_name", "qval"))
bulk_egenes <- bulk_table %>%
  filter(qval <= 0.05) %>%
  filter(gene_name %in% real$gene_name) %>%
  pull(gene_name)

# one row per gene: is it an eGene of each cell type, and of bulk? (TRUE / FALSE columns, as the UpSet needs)
all_genes <- unique(c(as.character(hits$gene_name), bulk_egenes))
membership <- data.frame(gene = all_genes)
for (ct in CELLTYPES) {
  genes_of_ct <- hits %>% filter(celltype == ct) %>% pull(gene_name)
  membership[[ct]] <- membership$gene %in% genes_of_ct
}
membership[["bulk (GTEx)"]] <- membership$gene %in% bulk_egenes
set_names <- c("bulk (GTEx)", rev(CELLTYPES))                        # order of the rows in the plot

p14 <- ComplexUpset::upset(
  membership,
  intersect = set_names,
  name = NULL,
  n_intersections = 30,
  sort_sets = FALSE,
  base_annotations = list(`genes in the intersection` = ComplexUpset::intersection_size(text = list(size = 2.6))),
  set_sizes = ComplexUpset::upset_set_size() + ylab("eGenes")) +
  labs(title = "Our cell-type eGenes and GTEx bulk whole-blood eGenes (all shown)")
p14
# the numbers behind it: per cell type, how many of its eGenes are bulk eGenes
share_bulk_all_tested <- mean(unique(real$gene_name) %in% bulk_egenes)  # the baseline: all tested genes
overlap_counts <- hits %>%
  mutate(is_bulk = gene_name %in% bulk_egenes) %>%
  group_by(celltype) %>%
  summarise(egenes = n(),
            also_bulk = sum(is_bulk),
            share_bulk = mean(is_bulk),
            new_vs_bulk = sum(!is_bulk),
            .groups = "drop") %>%
  mutate(share_bulk_all_tested_genes = share_bulk_all_tested)
print(as.data.frame(overlap_counts), digits = 2)
save_fig(p14, overlap_counts, "upsetplot_with_bulk", width = 14, height = 7)

## ================= 15. per-variant plots (step 10) of the example list =================
#this lets us view specific variants in specific cell types :)

# expression vs the cell type's proportion, one line per genotype, for every gene x cell type x variant in the list
# output/figures<SCAN_TAG>/variants/example_variants.csv (columns gene_name, celltype[, variant_id]; made by step 10b,
# or edit / replace it by hand to choose other hits). The plot function comes from 10_plot_variant.R, loaded into its
# own environment so it does not overwrite this script's objects.
example_file <- file.path(fig_dir, "variants", "example_variants.csv")
if (file.exists(example_file)) {
  
  # the list of examples (no variant given -> the gene's lead variant)
  examples <- read.csv(example_file)
  if (!"variant_id" %in% names(examples)) {
    examples$variant_id <- NA
  }
  
  # load plot_variant() from step 10 without running its own plots
  Sys.setenv(PLOT_VARIANT_NO_MAIN = "1")
  pv <- new.env()
  sys.source(file.path(pipe_dir, "10_plot_variant.R"), envir = pv)
  Sys.unsetenv("PLOT_VARIANT_NO_MAIN")
  
  # one plot per row; a row that cannot be plotted is skipped with a message
  for (r in seq_len(nrow(examples))) {
    gene     <- examples$gene_name[r]
    celltype <- examples$celltype[r]
    variant  <- examples$variant_id[r]
    p <- tryCatch(pv$plot_variant(gene, celltype, variant),
                  error = function(e) NULL)
    if (is.null(p)) {
      message("skipped ", gene, " in ", celltype)
    } else {
      all_plots[[paste0("09_15_variant_", gene, "_", celltype)]] <- p
    }
  }
  cat("per-variant plots:", nrow(examples), "rows of", example_file, "\n")
}

## ================= all figures in one PDF =================
#if we want to save all figures in one pdf later we can use this:
#pdf_file <- file.path(fig_dir, "09_all_figures.pdf")
#pdf(pdf_file, width = 13, height = 10)
#for (p in all_plots) print(p)
#invisible(dev.off())
#cat("all figures:", pdf_file, "\n")
#message("09_results_figures done")
