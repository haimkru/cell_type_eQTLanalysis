#author: haim krupkin+claude
#date: 11.07.2026 (11.08: also the raw-mixture run, for the robustness check)
#description: step 5. LM22's 22 cell types -> 6 broad classes, the proportions used by the eQTL step.
#   Neutrophils = Neutrophils
#   Monocytes   = Monocytes + Macrophages M0, M1, M2
#   Bcells      = B naive + B memory + Plasma
#   CD4         = CD4 naive, memory resting, memory activated, T follicular helper, Tregs
#   CD8         = CD8
#   NK          = NK resting + NK activated
# Left out (small): gamma delta T, dendritic cells, mast cells, eosinophils. The 6 classes are rescaled to sum to 1,
# with a floor of 1e-4 (the TCA model needs every proportion > 0).
# Done for both step-4 runs: main (batch-removed mixture) and raw (robustness check).
# Also prints the mean per class, pre- vs post-mortem (the batch removal should make these similar).
# Inputs : output/04_cibersortx/fractions.tsv, output/04_cibersortx_raw/fractions.tsv, output/02_covariates.csv
# Outputs: output/05_proportions_6classes.csv (main), output/05_proportions_6classes_raw.csv (check),
#          output/05_proportions_by_draw.csv
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
classes  <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
cov_df   <- read.csv(file.path(out_dir, "02_covariates.csv"))

## ================= 22 LM22 types -> 6 classes, floor, sum to 1 =================
to_6_classes <- function(fractions_file) {
  fr <- as.data.frame(fread(fractions_file))
  W <- data.frame(SAMPID = fr$Mixture)
  W$Neutrophils <- fr$Neutrophils
  W$Monocytes   <- fr$Monocytes + fr$Macrophages_M0 + fr$Macrophages_M1 + fr$Macrophages_M2
  W$Bcells      <- fr$B_cells_naive + fr$B_cells_memory + fr$Plasma_cells
  W$CD4         <- fr$T_cells_CD4_naive + fr$T_cells_CD4_memory_resting + fr$T_cells_CD4_memory_activated +
                   fr$T_cells_follicular_helper + fr$T_cells_regulatory_Tregs
  W$CD8         <- fr$T_cells_CD8
  W$NK          <- fr$NK_cells_resting + fr$NK_cells_activated
  m <- as.matrix(W[, classes])
  m <- pmax(m, 1e-4)
  m <- m / rowSums(m)
  W[, classes] <- m
  W
}
W_main <- to_6_classes(file.path(out_dir, "04_cibersortx/fractions.tsv"))
W_raw  <- to_6_classes(file.path(out_dir, "04_cibersortx_raw/fractions.tsv"))
round(colMeans(W_main[, classes]), 3)
round(colMeans(W_raw[, classes]), 3)

## ================= check: mean per class, pre- vs post-mortem, both runs =================
mean_by_draw <- function(W, label) {
  W %>%
    left_join(cov_df %>% select(SAMPID, post_mortem), by = "SAMPID") %>%
    mutate(draw = ifelse(post_mortem > 0.5, "post", "pre")) %>%
    pivot_longer(all_of(classes), names_to = "celltype", values_to = "fraction") %>%
    group_by(celltype, draw) %>%
    summarise(mean = mean(fraction), .groups = "drop") %>%
    pivot_wider(names_from = draw, values_from = mean) %>%
    mutate(run = label)
}
by_draw <- bind_rows(mean_by_draw(W_main, "main (batch removed)"), mean_by_draw(W_raw, "raw (check)"))
print(by_draw %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())

## ================= save =================
write.csv(W_main, file.path(out_dir, "05_proportions_6classes.csv"), row.names = FALSE)
write.csv(W_raw, file.path(out_dir, "05_proportions_6classes_raw.csv"), row.names = FALSE)
write.csv(by_draw, file.path(out_dir, "05_proportions_by_draw.csv"), row.names = FALSE)
out_files <- file.path(out_dir, c("05_proportions_6classes.csv", "05_proportions_6classes_raw.csv", "05_proportions_by_draw.csv"))
out_files
stopifnot(all(file.exists(out_files)))
message("05_proportions_6classes done")
