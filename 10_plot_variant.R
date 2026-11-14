#author: haim krupkin+claude
#date: 11.12.2026
#description: step 10. Plot any gene x variant x cell type, in the style of Kim-Hellmuth et al. 2020 (Science) Fig. 1C:
# expression (y) against the cell type's proportion (x), one colour and one regression line per genotype.
# A real cell-type eQTL shows lines that FAN OUT as the cell type grows (no genotype difference where the cell
# type makes none of the RNA, the largest where it makes most of it). Use it to pick good example variants.
#   y = log2(TPM + 1) with the known covariates regressed out (way of death, ischemic time, RIN, centre, sex, age,
#       5 genotype PCs; NOT the proportions, they are the x axis), then scaled to SD 1
#   x = the cell type's proportion (default); XAXIS=phi -> its share of the gene's RNA (step 5b; ~1 for marker genes,
#       so only useful for genes several cell types express)
# USAGE (from this folder):
#   one plot      : Rscript 10_plot_variant.R GENE CELLTYPE [VARIANT]      (no VARIANT: the gene's lead variant)
#                   e.g. Rscript 10_plot_variant.R KIR2DS4 NK
#   a list        : Rscript 10_plot_variant.R --file my_list.csv           (columns gene_name, celltype[, variant_id])
#   top hits      : Rscript 10_plot_variant.R --top 3                       (the 3 best eGenes of every cell type)
#   table to pick : Rscript 10_plot_variant.R --fan-table   -> figures<SCAN_TAG>/variants/10_fan_table.csv: every eGene
#                   with its fan (main effect, genotype x proportion interaction: beta, p) and the number of cell
#                   types that express the gene. Pick from it, then plot.
#   in RStudio    : set GENE / CELLTYPE / VARIANT in "what to plot" below, then run the whole script (Source, or
#                   ctrl+shift+enter): the plot shows in the Plots pane and is saved. Afterwards, in the console:
#                   MAIN FUNCTION: plot_gene() (all options in its comment block below)
#                     plot_gene("DCBLD1")                                NNLS expression + every cell type above the
#                                                                        RNA share threshold, same variant
#                     plot_gene("DCBLD1", celltypes = "all")             all 6 cell types
#                     plot_gene("DCBLD1", celltypes = c("NK"), profile = FALSE, scale = "linear", binned = TRUE)
#                     plot_variant("TRBV28", "CD8")                      a single cell type only
#                     plot_expression_profile("SERPINB10")               only the NNLS expression bars
#                     fan <- fan_table()                                 the table (takes a while), or read the CSV
#                     fan %>% filter(n_celltypes_expressing >= 2, p_int_prop < 1e-5) %>% arrange(desc(abs(fan_prop)))
# Style (11.13, set once in "plot style"): large text, no titles, little white space; legend (bottom) = genotypes,
# titled with the variant; on the plot: b and p of the genotype main effect and of the genotype x proportion interaction
# (simple fit expression ~ dosage * x, descriptive; the step-6 result is printed to the console).
# Inputs : output/01_tpm_noHb.rds, 02_covariates.csv, 05_proportions_6classes.csv, 05b_theta.csv,
#          output/genotypes/<chr>.{bed,bim,fam}, 06_scan_<chr>_top_variants.csv, figures/09_3_all_hits.csv
# Outputs: output/figures/variants/<gene>_<celltype>_<variant>.png (+ .csv of the plotted data);
#          for --file / --top also one PDF with all pages
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(snpStats))
suppressMessages(library(patchwork))     # plot_gene_panel: several plots in one figure

## ================= what to plot (RStudio) =================
GENE     <- "KIR2DS4"
CELLTYPE <- "NK"          # Neutrophils, Monocytes, Bcells, CD4, CD8, NK
VARIANT  <- NA            # NA = the gene's lead variant in this cell type; or e.g. "chr19_54815265_G_T_b38"

## ================= settings =================
pipe_dir <-"/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
plot_dir <- file.path(out_dir, paste0("figures", Sys.getenv("SCAN_TAG", "_epc20")), "variants")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
pc_file  <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
COVARS    <- c("post_mortem", "ventilator", "SMTSISCH", "SMRIN", "centre_B1", "male", "age")
XAXIS     <- Sys.getenv("XAXIS", "prop")                    # prop (cell proportion, default) or phi (RNA share)
SCAN_TAG  <- Sys.getenv("SCAN_TAG", "_epc20")               # 11.13: main result = 20 hidden-noise PCs; "" = the model without
GENO_COLS <- c("#1f77b4", "#ff7f0e", "#2ca02c")             # 0, 1, 2 copies of the ALT allele

## ================= plot style: large text, no titles, little white space (change it here for all plots) =================
BASE_SIZE  <- 20            # text size of axes and legends
LABEL_SIZE <- 6             # text written inside a plot (mm)
PANEL_SIZE <- 5             # inches per panel when saving (small file, large text): width
PANEL_HEIGHT <- 7           # inches per row: tall enough for the panel, the long axis titles and the legend below
DPI        <- 300
theme_big <- theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-8, 0, 0, 0),
        plot.margin = margin(4, 6, 4, 4),
        panel.grid.minor = element_blank())

## ================= data shared by all plots (read once) =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
W <- as.matrix(W_df[match(cov_df$SAMPID, W_df$SAMPID), CELLTYPES])
theta_df <- read.csv(file.path(out_dir, "05b_theta.csv"))
theta <- as.matrix(theta_df[, CELLTYPES])
rownames(theta) <- theta_df$gene_name
tpm <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
pcs <- fread(pc_file)
setnames(pcs, 1, "ID")
geno_pc <- t(as.matrix(pcs[ID %in% paste0("PC", 1:5), cov_df$SUBJID, with = FALSE]))
# the known covariates (not the proportions): regressed out of the expression
C <- cbind(1, scale(as.matrix(cov_df[, COVARS])), scale(geno_pc))
hits <- read.csv(file.path(out_dir, paste0("figures", SCAN_TAG), "09_3_all_hits.csv"))

# each gene's mean RNA share per cell type (step 5b): a cell type "expresses" the gene if its share >= 0.10
# (the same gate as step 6: only those cell types were tested for this gene)
PHI_GATE <- 0.10
mean_phi <- read.csv(file.path(out_dir, "05b_mean_phi.csv"))

## ================= helpers =================
# caches: each chromosome's variant list and scan results are read once, not once per plot
bim_cache <- list()
top_cache <- list()

get_bim <- function(chr) {
  if (is.null(bim_cache[[chr]])) {
    bim_cache[[chr]] <<- fread(file.path(out_dir, "genotypes", paste0(chr, ".bim")), header = FALSE)
  }
  bim_cache[[chr]]
}

get_top <- function(chr) {
  if (is.null(top_cache[[chr]])) {
    top <- fread(file.path(out_dir, paste0("06_scan_", chr, "_top_variants", SCAN_TAG, ".csv")))
    top_cache[[chr]] <<- top[set == "genotype"]
  }
  top_cache[[chr]]
}

# chromosome of a gene: from the hit list, else from the variant ID, else search the step-6 scans
find_chr <- function(gene, variant = NA) {
  if (gene %in% hits$gene_name) return(hits$chr[hits$gene_name == gene][1])
  if (!is.na(variant)) return(sub("_.*", "", variant))
  for (ch in paste0("chr", 1:22)) {
    g <- fread(file.path(out_dir, paste0("06_scan_", ch, "_genes", SCAN_TAG, ".csv")), select = "gene_name")
    if (gene %in% g$gene_name) return(ch)
  }
  stop("gene not found in the step-6 scans: ", gene)
}

# the scan's row for this gene x cell type (real genotypes): lead variant, effects, p-values
scan_row <- function(gene, celltype, chr) {
  top <- get_top(chr)
  top[gene_name == gene & test == celltype]
}

# dosages read in advance (prefetch_dosages): reading many variants of a chromosome at once is much faster than one by
# one (read.plink re-reads the whole variant list on every call)
dosage_cache <- list()

prefetch_dosages <- function(variants, chr) {
  pre <- file.path(out_dir, "genotypes", chr)
  bim <- get_bim(chr)
  k <- which(bim$V2 %in% unique(variants))
  if (length(k) == 0) return(invisible(NULL))
  pl <- read.plink(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"), select.snps = k)
  G <- as(pl$genotypes, "numeric")
  alt <- sapply(strsplit(colnames(G), "_"), `[`, 4)
  flip <- pl$map$allele.2 != alt
  G[, flip] <- 2 - G[, flip]                               # count of the ALT allele
  G <- G[match(cov_df$SUBJID, rownames(G)), , drop = FALSE]
  for (v in colnames(G)) {
    dosage_cache[[v]] <<- G[, v]
  }
  invisible(NULL)
}

# ALT-allele dosage of one variant for the 670 donors (same coding as step 6)
read_dosage <- function(variant, chr) {
  if (!is.null(dosage_cache[[variant]])) return(dosage_cache[[variant]])
  pre <- file.path(out_dir, "genotypes", chr)
  bim <- get_bim(chr)
  k <- which(bim$V2 == variant)
  if (length(k) != 1) stop("variant not in the genotype file: ", variant)
  pl <- read.plink(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"), select.snps = k)
  g <- as(pl$genotypes, "numeric")[, 1]
  alt <- strsplit(variant, "_")[[1]][4]
  if (pl$map$allele.2 != alt) g <- 2 - g                   # count of the ALT allele
  g[match(cov_df$SUBJID, names(g))]
}

# in which cell types is the gene expressed (mean RNA share >= PHI_GATE)?
celltypes_expressing <- function(gene) {
  row <- mean_phi %>% filter(gene_name == gene)
  if (nrow(row) == 0) return(character(0))
  shares <- unlist(row[1, CELLTYPES])
  CELLTYPES[shares >= PHI_GATE]
}

## ================= the data behind one plot =================
# one row per donor: proportion, RNA share (phi), expression (covariates removed, SD units), dosage, genotype
# scale = "log"   : log2(TPM + 1), covariates removed, scaled to SD units (the scale of the step-6 model)
# scale = "linear": TPM, covariates removed (residuals + the gene's mean TPM). On this scale an effect in one cell
#                   type grows with that cell type's amount (the lines fan out from ~0); on the log scale the same
#                   effect is a shift when the cell type makes most of the gene's RNA.
variant_data <- function(gene, celltype, variant, chr, scale = "log") {
  stopifnot(scale %in% c("log", "linear"))
  tpm_gene <- as.numeric(tpm[gene, cov_df$SAMPID])
  if (scale == "log") {
    y <- log2(tpm_gene + 1)
    y <- as.numeric(scale(lm.fit(C, y)$residuals))
  } else {
    y <- lm.fit(C, tpm_gene)$residuals + mean(tpm_gene)
  }
  # the cell type's share of this gene's RNA per donor (phi), from the proportions and the gene's theta
  num <- sweep(W, 2, theta[gene, ], "*")
  phi <- (num / rowSums(num))[, celltype]
  # genotype labels from the variant ID (chr_pos_REF_ALT_b38)
  g <- read_dosage(variant, chr)
  ref <- strsplit(variant, "_")[[1]][3]
  alt <- strsplit(variant, "_")[[1]][4]
  labels <- c(paste0(ref, ref), paste0(ref, alt), paste0(alt, alt))
  d <- data.frame(SAMPID = cov_df$SAMPID, prop = W[, celltype], phi = phi, expression = y, dosage = g)
  d <- d %>% filter(!is.na(dosage))
  d <- d %>% mutate(genotype = factor(labels[round(dosage) + 1], levels = labels))
  d
}

## ================= main effect, interaction, and how much the lines fan out =================
# fit: expression ~ dosage * x  (x = the proportion or phi)
#   beta_main / p_main : the genotype (dosage) term = per-allele effect where x = 0
#   beta_int  / p_int  : the dosage:x term = how the per-allele effect changes per unit of x (the fan)
#   effect_low / effect_high: the per-allele effect (SD units) at the 10th / 90th percentile of x
#   fan = effect_high - effect_low  (0 = parallel lines; large = the lines spread as the cell type grows)
# For a gene made by ONE cell type only, phi is ~1 for everyone (no variation in x: NA), and on the proportion axis the
# lines are about parallel (a shift, not a fan): its cell-type specificity is in the expression, not the regulation.
# NOTE: a simple descriptive fit (no other proportions, no hidden-noise PCs); the eQTL test is step 6's joint model.
fan_stats <- function(d, x_column) {
  empty <- data.frame(beta_main = NA, p_main = NA, beta_int = NA, p_int = NA,
                      effect_low = NA, effect_high = NA, fan = NA)
  x <- d[[x_column]]
  if (sd(x) < 1e-3) return(empty)                       # x hardly varies: no fan to measure
  fit_data <- data.frame(expression = d$expression, dosage = d$dosage, x = x)
  fit <- lm(expression ~ dosage * x, data = fit_data)
  coefs <- summary(fit)$coefficients
  if (!"dosage:x" %in% rownames(coefs)) return(empty)
  x_low  <- unname(quantile(x, 0.1))
  x_high <- unname(quantile(x, 0.9))
  effect_low  <- coefs["dosage", "Estimate"] + coefs["dosage:x", "Estimate"] * x_low
  effect_high <- coefs["dosage", "Estimate"] + coefs["dosage:x", "Estimate"] * x_high
  data.frame(beta_main = coefs["dosage", "Estimate"],
             p_main    = coefs["dosage", "Pr(>|t|)"],
             beta_int  = coefs["dosage:x", "Estimate"],
             p_int     = coefs["dosage:x", "Pr(>|t|)"],
             effect_low = effect_low,
             effect_high = effect_high,
             fan = effect_high - effect_low)
}

## ================= the fan as one line: the per-allele effect in bins of x =================
# donors split into n_bins equal groups by x; in each group: expression ~ dosage -> per-allele effect and its SE.
# A cell-type effect: ~0 in the lowest group, rising (or falling) steadily with x.
binned_effect <- function(d, n_bins = 5) {
  d$bin <- ntile(d$x, n_bins)
  rows <- list()
  for (k in seq_len(n_bins)) {
    d_bin <- d %>% filter(bin == k)
    coefs <- summary(lm(expression ~ dosage, data = d_bin))$coefficients
    rows[[k]] <- data.frame(bin = k,
                            x_median = median(d_bin$x),
                            donors = nrow(d_bin),
                            effect = coefs["dosage", "Estimate"],
                            se = coefs["dosage", "Std. Error"])
  }
  bind_rows(rows)
}

plot_binned_effect <- function(d, x_label, y_label, n_bins = 5) {
  b <- binned_effect(d, n_bins)
  ggplot(b, aes(x_median, effect)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_line(colour = "grey30", linewidth = 1) +
    geom_pointrange(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), colour = "black",
                    size = 0.8, linewidth = 1) +
    labs(x = x_label, y = y_label) +   # one point per fifth of donors: mean +- 95% CI
    theme_big
}

## ================= one plot =================
# xaxis : "prop" (the cell type's proportion) or "phi" (its share of the gene's RNA)
# scale : "log" (log2 TPM, SD units; the model's scale) or "linear" (TPM; shows the fan)
# binned: TRUE adds the per-allele effect in fifths of x next to the scatter
plot_variant <- function(gene, celltype, variant = NA, save = TRUE, xaxis = XAXIS, scale = "log", binned = FALSE,
                         parts = FALSE) {
  stopifnot(celltype %in% CELLTYPES, gene %in% rownames(tpm), gene %in% rownames(theta), xaxis %in% c("prop", "phi"))
  chr <- find_chr(gene, variant)
  sr <- scan_row(gene, celltype, chr)
  if (is.na(variant)) {
    if (nrow(sr) == 0) stop("no scan result for ", gene, " in ", celltype, " (cell type gated out?): give a VARIANT")
    variant <- sr$variant_id
  }
  # the data and the main / interaction fit
  d <- variant_data(gene, celltype, variant, chr, scale = scale)
  d$x <- d[[xaxis]]
  st <- fan_stats(d, xaxis)
  expressing <- celltypes_expressing(gene)

  # text written on the plot: main effect and interaction (beta, p)
  # main = genotype term; interaction = genotype x x-axis term (short, so it fits in a small panel)
  stats_text <- paste0(sprintf("main: b = %.2f, p = %.1e", st$beta_main, st$p_main), "\n",
                       sprintf("interaction: b = %.2f, p = %.1e", st$beta_int, st$p_int))

  # printed to the console (not on the plot): the step-6 result at this variant, expressing cell types, n per genotype
  # (any saved lead-variant row of this gene at this variant carries all 6 effects)
  rows_at_variant <- get_top(chr)[gene_name == gene & variant_id == variant]
  has_scan <- nrow(rows_at_variant) >= 1 && !is.na(rows_at_variant[[paste0("b_", celltype)]][1])
  if (has_scan) {
    cat(sprintf("%s %s %s: joint model (step 6) effect %.2f (SE %.2f), p = %.2g\n", gene, celltype, variant,
                rows_at_variant[[paste0("b_", celltype)]][1], rows_at_variant[[paste0("se_", celltype)]][1],
                rows_at_variant[[paste0("p_", celltype)]][1]))
  }
  cat("  expressed in:", paste(expressing, collapse = ", "), "; n per genotype:",
      paste(table(d$genotype), collapse = " / "), "\n")

  labels <- levels(d$genotype)
  x_label <- if (xaxis == "prop") paste(celltype, "proportion") else paste(celltype, "RNA share")
  y_label <- if (scale == "log") paste(gene, "(log2 TPM, SD)") else paste(gene, "(TPM)")
  p <- ggplot(d, aes(x, expression, colour = genotype, fill = genotype)) +
    geom_point(size = 1.5, alpha = 0.6) +
    # one line per genotype with its 95% confidence band (shaded)
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95, linewidth = 1.8, alpha = 0.2) +
    annotate("label", x = -Inf, y = Inf, label = stats_text, hjust = -0.02, vjust = 1.05, size = LABEL_SIZE * 0.65,
             fill = "white") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +     # room at the top for the label
    scale_colour_manual(values = setNames(GENO_COLS, labels), drop = FALSE) +
    scale_fill_manual(values = setNames(GENO_COLS, labels), drop = FALSE) +
    labs(x = x_label, y = y_label, colour = variant, fill = variant) +   # legend title = the variant (one legend)
    guides(colour = guide_legend(title.position = "top", override.aes = list(size = 4)),
           fill = guide_legend(title.position = "top")) +   # same title + labels -> one legend
    theme_big
  # optional: the binned per-allele effect next to the scatter
  width <- PANEL_SIZE
  p_binned <- NULL
  if (binned) {
    effect_label <- if (scale == "log") "per-allele effect (SD)" else "per-allele effect (TPM)"
    p_binned <- plot_binned_effect(d, x_label, effect_label)
  }
  # parts = TRUE: return the scatter and the binned plot separately (plot_gene arranges them)
  if (parts) return(list(scatter = p, binned = p_binned))
  if (binned) {
    p <- p + p_binned + plot_layout(ncol = 2)
    width <- 2 * PANEL_SIZE
  }
  if (save) {
    base <- file.path(plot_dir, paste0(gene, "_", celltype, "_", variant, if (xaxis == "phi") "_phi",
                                       if (scale == "linear") "_linear", if (binned) "_binned"))
    ggsave(paste0(base, ".png"), p, width = width, height = PANEL_HEIGHT, dpi = DPI)
    write.csv(d, paste0(base, ".csv"), row.names = FALSE)
    cat("figure:", paste0(base, ".png"), "\n")
  }
  p
}

## ================= what the joint model sees: partial residuals =================
# The step-6 model refitted for ONE gene x variant, exactly as step 6 does it:
#   y = log2(TPM + 1); weighted by 1 / (expected SD of each sample, from the step-6 null fit);
#   design = W | W x death | covariates (RIN, ischemic time, centre, sex, age, 5 genotype PCs, 20 hidden-noise PCs)
#            | genotype x phi of each cell type above the RNA share threshold | genotype x death;
#   robust (HC3) standard errors, p from the t distribution -> the same p as step 6 (checked, printed).
# Partial residual for cell type c = what is left after removing everything the model fitted EXCEPT c's genotype term:
#   residual + b_c x phi_c x g. Plotted against phi_c, the genotype lines start together at phi = 0 and spread with phi:
#   this is the signal the test for c uses, and nothing else (no other cell types, no covariates, no noise PCs).
EPC_FILE <- file.path(out_dir, "05c_expression_pcs.csv")
N_EPC <- 20
OTHER_COV <- c("SMTSISCH", "SMRIN", "centre_B1", "male", "age")
DEATH <- c("post_mortem", "ventilator")
joint_cache <- list()

# the base design of step 6 (same for every gene); built once
get_joint_base <- function() {
  if (!is.null(joint_cache$base)) return(joint_cache)
  death <- scale(as.matrix(cov_df[, DEATH]))
  epc_df <- read.csv(EPC_FILE)
  epc <- as.matrix(epc_df[match(cov_df$SAMPID, epc_df$SAMPID), paste0("ePC", 1:N_EPC)])
  C2 <- scale(cbind(as.matrix(cov_df[, OTHER_COV]), geno_pc))
  C2 <- scale(cbind(C2, epc))
  W_death <- cbind(W * death[, 1], W * death[, 2])
  colnames(W_death) <- c(paste0(CELLTYPES, "_x_", DEATH[1]), paste0(CELLTYPES, "_x_", DEATH[2]))
  colnames(C2) <- paste0("cov", seq_len(ncol(C2)))
  joint_cache$death <<- death
  joint_cache$base  <<- cbind(W, W_death, C2)
  joint_cache
}

# the step-6 null fit of this chromosome: each gene's SD per cell type (sigmas) and tau -> the sample weights
get_null_fit <- function(chr) {
  key <- paste0("null_", chr)
  if (is.null(joint_cache[[key]])) {
    joint_cache[[key]] <<- readRDS(file.path(out_dir, paste0("06_scan_", chr, SCAN_TAG, "_null_fit.rds")))
  }
  joint_cache[[key]]
}

# plant_celltype / plant_sd (optional, for the spike-in simulation, as step 6 MODEL=spike): before fitting, add a known
# effect in ONE cell type: y + phi_c x g x plant_sd x SD(y). The true effect (log2 units) is returned as true_effect.
joint_refit <- function(gene, variant, chr, plant_celltype = NA, plant_sd = 0) {
  jc <- get_joint_base()
  # sample weights: 1 / expected SD of each sample
  nf <- get_null_fit(chr)
  sd_i <- sqrt(drop(W^2 %*% nf$sigmas[gene, ]^2) + nf$tau^2)
  # outcome and RNA shares
  y <- log2(as.numeric(tpm[gene, cov_df$SAMPID]) + 1)
  num <- sweep(W, 2, theta[gene, ], "*")
  phi <- num / rowSums(num)
  active <- colMeans(phi) >= PHI_GATE
  # centred genotype (missing -> mean, as step 6)
  g <- read_dosage(variant, chr)
  g[is.na(g)] <- mean(g, na.rm = TRUE)
  g_c <- g - mean(g)
  # spike-in: plant a known effect in one cell type (through its RNA share, exactly as step 6 does)
  true_effect <- 0
  if (!is.na(plant_celltype)) {
    true_effect <- plant_sd * sd(y)
    y <- y + phi[, plant_celltype] * g_c * true_effect
  }
  # genotype columns: x phi of each active cell type, x each death variable
  geno_cols <- cbind(phi[, active, drop = FALSE] * g_c, jc$death * g_c)
  colnames(geno_cols) <- c(paste0("g_", CELLTYPES[active]), paste0("g_x_", DEATH))
  X_all <- cbind(jc$base, geno_cols)
  # weighted least squares = ordinary least squares on the scaled rows
  Xs <- X_all / sd_i
  ys <- y / sd_i
  XtX_inv <- solve(crossprod(Xs))
  b <- drop(XtX_inv %*% crossprod(Xs, ys))
  e <- drop(ys - Xs %*% b)
  # HC3 robust covariance: each sample's squared residual inflated by its leverage
  h <- rowSums((Xs %*% XtX_inv) * Xs)
  meat <- crossprod(Xs * (e / (1 - h)))
  V <- XtX_inv %*% meat %*% XtX_inv
  se <- sqrt(diag(V))
  df <- length(y) - ncol(Xs)
  p <- 2 * pt(-abs(b / se), df)
  names(b) <- names(se) <- names(p) <- colnames(X_all)
  list(b = b, se = se, p = p, resid = e * sd_i, phi = phi, g = g, g_c = g_c, active = active,
       true_effect = true_effect)
}

plot_partial_residual <- function(gene, celltype, variant = NA, save = TRUE) {
  chr <- find_chr(gene, variant)
  if (is.na(variant)) variant <- scan_row(gene, celltype, chr)$variant_id
  fit <- joint_refit(gene, variant, chr)
  term <- paste0("g_", celltype)
  if (!term %in% names(fit$b)) stop(celltype, " is under the RNA share threshold for ", gene, ": not in the model")
  b_c <- fit$b[[term]]
  se_c <- fit$se[[term]]
  p_c <- fit$p[[term]]
  # check: the same as step 6?
  rows_at_variant <- get_top(chr)[gene_name == gene & variant_id == variant]
  if (nrow(rows_at_variant) >= 1) {
    cat(sprintf("%s %s: refit b = %.3f, SE = %.3f, p = %.2g | step 6: b = %.3f, SE = %.3f, p = %.2g\n", gene, celltype,
                b_c, se_c, p_c, rows_at_variant[[paste0("b_", celltype)]][1],
                rows_at_variant[[paste0("se_", celltype)]][1], rows_at_variant[[paste0("p_", celltype)]][1]))
  }
  # partial residual: residual + this cell type's genotype term
  ref <- strsplit(variant, "_")[[1]][3]
  alt <- strsplit(variant, "_")[[1]][4]
  labels <- c(paste0(ref, ref), paste0(ref, alt), paste0(alt, alt))
  # x: the cell type's RNA share (what the model uses). A gene made by this cell type only has share ~1 in everyone
  # (nothing to spread over), so then the cell PROPORTION is used instead.
  share <- fit$phi[, celltype]
  use_prop <- sd(share) < 0.02
  x_values <- if (use_prop) W[, celltype] else share
  x_label  <- if (use_prop) paste(celltype, "proportion") else paste(celltype, "RNA share")
  d <- data.frame(x = x_values,
                  partial = fit$resid + b_c * share * fit$g_c,
                  genotype = factor(labels[round(fit$g) + 1], levels = labels))
  stats_text <- sprintf("joint model: b = %.2f, p = %.1e", b_c, p_c)
  p <- ggplot(d, aes(x, partial, colour = genotype, fill = genotype)) +
    geom_point(size = 1.5, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1.8, alpha = 0.2) +
    annotate("label", x = -Inf, y = Inf, label = stats_text, hjust = -0.02, vjust = 1.1, size = LABEL_SIZE * 0.65,
             fill = "white") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    scale_colour_manual(values = setNames(GENO_COLS, labels), drop = FALSE) +
    scale_fill_manual(values = setNames(GENO_COLS, labels), drop = FALSE) +
    labs(x = x_label, y = paste(gene, "partial resid."),
         colour = variant, fill = variant) +
    guides(colour = guide_legend(title.position = "top"), fill = guide_legend(title.position = "top")) +
    theme_big
  if (save) {
    out_file <- file.path(plot_dir, paste0("partial_", gene, "_", celltype, "_", variant, ".png"))
    ggsave(out_file, p, width = PANEL_SIZE, height = PANEL_HEIGHT, dpi = DPI)
    cat("figure:", out_file, "\n")
  }
  p
}

## ================= expression of the gene in each cell type (NNLS, step 5b) =================
# bars: theta = the gene's estimated expression per cell type (TPM, from regressing bulk TPM on the 6 proportions;
# 0.001 = not expressed). Label: the cell type's mean share of the gene's RNA. Dark = expressed (share >= PHI_GATE:
# the cell types step 6 tested for this gene).
plot_expression_profile <- function(gene) {
  stopifnot(gene %in% rownames(theta))
  share_row <- mean_phi %>% filter(gene_name == gene)
  shares <- if (nrow(share_row) == 1) as.numeric(unlist(share_row[1, CELLTYPES])) else rep(NA, length(CELLTYPES))
  prof <- data.frame(celltype = factor(CELLTYPES, levels = CELLTYPES),
                     theta = as.numeric(theta[gene, CELLTYPES]),
                     share = shares)
  prof <- prof %>%
    mutate(status = ifelse(share >= PHI_GATE, "Above RNA share threshold", "Under RNA share threshold"),
           label = paste0(round(100 * share), "%"))       # bar label: mean share of the gene's RNA
  ggplot(prof, aes(celltype, theta, fill = status)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = label), vjust = -0.3, size = LABEL_SIZE) +
    scale_fill_manual(values = c("Above RNA share threshold" = "grey20", "Under RNA share threshold" = "grey75")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = paste(gene, "TPM (NNLS)"), fill = NULL) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_big +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

## ================= which cell type carries the effect? joint-model effects (+ OneK1K) =================
# The scatter plots fit one cell type at a time, so a cell type that shares the gene with the real one shows the effect
# too. The joint model separates them: this plot shows its 6 effects at one variant (+- 95% CI); cell types under the
# RNA share threshold were not tested. Optional second panel: OneK1K (single cells, other donors) - for each OneK1K
# cell type, -log10 q of the gene's top eQTL (each type's own lead variant; "not tested" = gene not expressed there).
ONEK1K_DIR <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/onek1k_top_eqtl/unzipped/OneK1K_TensorQTL_top_eQTL_summary"
ONEK1K_TYPES <- c("Mono_C", "Mono_NC", "DC", "B_IN", "B_MEM", "Plasma", "CD4_NC", "CD4_ET", "CD4_SOX4",
                  "CD8_NC", "CD8_ET", "CD8_S100B", "NK", "NK_R")

plot_joint_effects <- function(gene, variant, chr) {
  row <- get_top(chr)[gene_name == gene & variant_id == variant]
  if (nrow(row) == 0) stop("no step-6 result for ", gene, " at ", variant)
  # one row per cell type: effect, SE, p (NA = not tested: under the RNA share threshold)
  eff <- data.frame(celltype = CELLTYPES,
                    effect = as.numeric(unlist(row[1, paste0("b_", CELLTYPES), with = FALSE])),
                    se     = as.numeric(unlist(row[1, paste0("se_", CELLTYPES), with = FALSE])),
                    p      = as.numeric(unlist(row[1, paste0("p_", CELLTYPES), with = FALSE])))
  eff <- eff %>%
    mutate(celltype = factor(celltype, levels = rev(CELLTYPES)),
           status = case_when(is.na(effect) ~ "not tested",
                              p < 0.05 ~ "p < 0.05",
                              TRUE ~ "p >= 0.05"),
           p_label = ifelse(is.na(p), "not tested", paste0("p = ", signif(p, 2))))
  ggplot(eff, aes(effect, celltype, colour = status)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = effect - 1.96 * se, xmax = effect + 1.96 * se), width = 0.25, linewidth = 1.2,
                  orientation = "y", na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    geom_text(aes(x = -Inf, label = p_label), data = eff %>% filter(status == "not tested"),
              hjust = -0.1, colour = "grey60", size = LABEL_SIZE * 0.7) +     # at the left edge, always fully visible
    geom_text(aes(label = p_label), data = eff %>% filter(status != "not tested"),
              vjust = -1.1, size = LABEL_SIZE * 0.7, show.legend = FALSE) +
    scale_colour_manual(values = c("p < 0.05" = "black", "p >= 0.05" = "grey55", "not tested" = "grey80")) +
    labs(x = "joint-model effect (per allele, 95% CI)", y = NULL, colour = NULL) +
    theme_big
}

plot_onek1k_effects <- function(gene, chr) {
  rows <- list()
  for (ty in ONEK1K_TYPES) {
    f <- file.path(ONEK1K_DIR, paste0("OneK1K_", ty, ".sig_cis_qtl_pairs.", chr, ".csv"))
    q <- NA
    if (file.exists(f)) {
      x <- fread(f, select = c("phenotype_id", "qval"))
      x <- x[phenotype_id == gene]
      if (nrow(x) > 0) q <- min(x$qval)
    }
    rows[[ty]] <- data.frame(type = ty, q = q)
  }
  ok <- bind_rows(rows)
  # the gene is in no OneK1K cell type (e.g. antibody V genes are not in OneK1K's gene set): no panel
  if (all(is.na(ok$q))) {
    cat(gene, "is not in any OneK1K cell type: OneK1K panel left out\n")
    return(NULL)
  }
  ok <- ok %>%
    mutate(type = factor(type, levels = rev(ONEK1K_TYPES)),
           score = ifelse(is.na(q), 0, -log10(q)),
           status = case_when(is.na(q) ~ "not tested",
                              q < 0.05 ~ "q < 0.05",
                              TRUE ~ "q >= 0.05"))
  ggplot(ok, aes(score, type, fill = status)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    geom_text(aes(x = 0, label = "not tested"), data = ok %>% filter(status == "not tested"),
              hjust = -0.05, colour = "grey60", size = LABEL_SIZE * 0.6) +
    scale_fill_manual(values = c("q < 0.05" = "black", "q >= 0.05" = "grey65", "not tested" = "grey90")) +
    labs(x = "OneK1K eQTL, -log10 q", y = NULL, fill = NULL) +
    theme_big
}

# NNLS expression | joint-model effects | OneK1K (| optional: the scatter of one cell type)
# Example: plot_effects("KLRC3")  or  plot_effects("KLRC3", scatter_celltype = "NK")
#          plot_effects("KLRC3", partial_celltype = "NK")   last panel: partial residuals (p = the joint model's)
#          plot_effects("KLRC3", partial_celltype = "expressed")   partial residuals of every tested cell type (CD8 + NK)
plot_effects <- function(gene, variant = NA, onek1k = TRUE, scatter_celltype = NA, partial_celltype = NA,
                         save = TRUE) {
  chr <- find_chr(gene, variant)
  # the variant: the lead of the gene's strongest eGene call
  if (is.na(variant)) {
    gene_hits <- hits %>% filter(gene_name == gene) %>% arrange(p_gene)
    if (nrow(gene_hits) == 0) stop(gene, " is not an eGene: give a VARIANT")
    variant <- scan_row(gene, gene_hits$celltype[1], chr)$variant_id
  }
  plots <- list(plot_expression_profile(gene), plot_joint_effects(gene, variant, chr))
  if (onek1k) {
    p_onek1k <- plot_onek1k_effects(gene, chr)
    if (!is.null(p_onek1k)) plots <- c(plots, list(p_onek1k))
  }
  if (!is.na(scatter_celltype)) plots <- c(plots, list(plot_variant(gene, scatter_celltype, variant, save = FALSE)))
  # partial residuals: one panel per cell type; "expressed" = every cell type above the RNA share threshold
  if (identical(partial_celltype, "expressed")) partial_celltype <- celltypes_expressing(gene)
  if (!all(is.na(partial_celltype))) {
    for (ct in partial_celltype) {
      plots <- c(plots, list(plot_partial_residual(gene, ct, variant, save = FALSE)))
    }
  }
  combined <- wrap_plots(plots, nrow = 1)
  if (save) {
    out_file <- file.path(plot_dir, paste0("effects_", gene, "_", variant, ".png"))
    ggsave(out_file, combined, width = PANEL_SIZE * length(plots), height = PANEL_HEIGHT, dpi = DPI)
    cat("figure:", out_file, "  (variant:", variant, ")\n")
  }
  combined
}

## ================= THE plotting function: one gene, any cell types, optionally its expression profile =================
# gene      : the gene
# variant   : NA = the lead variant of the gene's strongest eGene call (the same variant in every panel)
# celltypes : "expressed" = the cell types above the RNA share threshold (PHI_GATE; the ones step 6 tested)
#             "all"       = all 6 cell types
#             or a vector, e.g. c("NK", "CD8")
# profile   : TRUE = first panel: the gene's NNLS expression per cell type (step 5b)
# xaxis     : "prop" (cell proportion) or "phi" (RNA share);  scale: "log" or "linear" (TPM; shows the fan)
# binned    : TRUE = each cell type also gets the per-allele effect in fifths of donors
# ncol      : panels per row (NA = up to 3)
# Examples: plot_gene("DCBLD1")                                   profile + the 3 expressing cell types
#           plot_gene("DCBLD1", celltypes = "all", profile = FALSE)  all 6 cell types
#           plot_gene("PADI2", celltypes = "Neutrophils", scale = "linear", binned = TRUE)
plot_gene <- function(gene, variant = NA, celltypes = "expressed", profile = TRUE, xaxis = XAXIS,
                      scale = "log", binned = FALSE, ncol = NA, save = TRUE) {
  chr <- find_chr(gene, variant)
  # the variant
  if (is.na(variant)) {
    gene_hits <- hits %>% filter(gene_name == gene) %>% arrange(p_gene)
    if (nrow(gene_hits) == 0) stop(gene, " is not an eGene: give a VARIANT")
    sr <- scan_row(gene, gene_hits$celltype[1], chr)
    variant <- sr$variant_id
  }
  # the cell types to show
  if (identical(celltypes, "expressed")) celltypes <- celltypes_expressing(gene)
  if (identical(celltypes, "all")) celltypes <- CELLTYPES
  if (length(celltypes) == 0) stop(gene, ": no cell type above the RNA share threshold (", PHI_GATE, ")")
  # the panels (all equal size)
  plots <- list()
  if (!binned) {
    # profile, then one scatter per cell type; up to 3 per row
    if (profile) plots[["profile"]] <- plot_expression_profile(gene)
    for (ct in celltypes) {
      plots[[ct]] <- plot_variant(gene, ct, variant, save = FALSE, xaxis = xaxis, scale = scale)
    }
    if (is.na(ncol)) ncol <- min(length(plots), 3)
  } else {
    # one row per cell type: [profile (first row only; empty below)] [scatter] [binned effect]
    for (i in seq_along(celltypes)) {
      ct <- celltypes[i]
      if (profile) {
        first_panel <- if (i == 1) plot_expression_profile(gene) else plot_spacer()
        plots[[paste0("first_", ct)]] <- first_panel
      }
      both <- plot_variant(gene, ct, variant, save = FALSE, xaxis = xaxis, scale = scale, binned = TRUE, parts = TRUE)
      plots[[paste0("scatter_", ct)]] <- both$scatter
      plots[[paste0("binned_", ct)]] <- both$binned
    }
    ncol <- if (profile) 3 else 2
  }
  combined <- wrap_plots(plots, ncol = ncol)
  if (save) {
    n_rows <- ceiling(length(plots) / ncol)
    out_file <- file.path(plot_dir, paste0("gene_", gene, "_", variant, if (xaxis == "phi") "_phi",
                                           if (scale == "linear") "_linear", if (binned) "_binned", ".png"))
    ggsave(out_file, combined, width = PANEL_SIZE * ncol, height = PANEL_HEIGHT * n_rows, dpi = DPI)
    cat("figure:", out_file, "\n")
  }
  combined
}
# old name, kept so earlier calls still work
plot_gene_panel <- function(gene, variant = NA, celltypes = NA, xaxis = XAXIS, save = TRUE, scale = "log") {
  if (all(is.na(celltypes))) celltypes <- "expressed"
  plot_gene(gene, variant, celltypes = celltypes, xaxis = xaxis, scale = scale, save = save)
}

## ================= the table to choose from: fan and expression for every eGene =================
# one row per eGene x cell type (09_3_all_hits.csv), at the scan's lead variant:
#   n_celltypes_expressing / celltypes_expressing : cell types with RNA share >= PHI_GATE (1 = made by one cell type)
#   z_joint       : the step-6 effect / SE in this cell type
#   *_prop, *_phi : main effect, interaction and fan with x = the proportion, and with x = the RNA share
#   min_genotype_n: donors in the smallest genotype group (small = a noisy line)
# Use it to choose, then plot_variant(gene, celltype). Saved to figures<SCAN_TAG>/variants/10_fan_table.csv
fan_table <- function(todo = hits, save = TRUE) {
  # read the genotypes of all the lead variants first, one chromosome at a time (much faster than one by one)
  for (ch in unique(todo$chr)) {
    genes_on_chr <- todo$gene_name[todo$chr == ch]
    lead_variants <- get_top(ch) %>%
      filter(gene_name %in% genes_on_chr) %>%
      pull(variant_id)
    prefetch_dosages(lead_variants, ch)
    cat("genotypes read:", ch, "\n")
  }
  rows <- list()
  for (r in seq_len(nrow(todo))) {
    gene     <- todo$gene_name[r]
    celltype <- todo$celltype[r]
    res <- tryCatch({
      chr <- find_chr(gene)
      sr  <- scan_row(gene, celltype, chr)
      if (nrow(sr) != 1) stop("no scan row")
      variant <- sr$variant_id
      d <- variant_data(gene, celltype, variant, chr)
      st_prop <- fan_stats(d, "prop")
      st_phi  <- fan_stats(d, "phi")
      names(st_prop) <- paste0(names(st_prop), "_prop")
      names(st_phi)  <- paste0(names(st_phi), "_phi")
      expressing <- celltypes_expressing(gene)
      info <- data.frame(gene_name = gene, celltype = celltype, variant_id = variant, chr = chr,
                         n_celltypes_expressing = length(expressing),
                         celltypes_expressing = paste(expressing, collapse = ","),
                         mean_phi = mean(d$phi),
                         z_joint = sr[[paste0("b_", celltype)]] / sr[[paste0("se_", celltype)]],
                         min_genotype_n = min(table(d$genotype)))
      cbind(info, st_prop, st_phi)
    }, error = function(e) {
      message("skipped ", gene, " ", celltype, ": ", conditionMessage(e))
      NULL
    })
    rows[[r]] <- res
    if (r %% 200 == 0) cat(r, "of", nrow(todo), "\n")
  }
  tab <- bind_rows(rows)
  # keep the lineage-conflict flag of the hit list, if it is there
  if ("lineage_conflict" %in% names(todo)) {
    flags <- todo %>% select(gene_name, celltype, lineage_conflict)
    tab <- tab %>% left_join(flags, by = c("gene_name", "celltype"))
  }
  tab <- tab %>% arrange(celltype, p_int_prop)
  if (save) {
    out_file <- file.path(plot_dir, "10_fan_table.csv")
    write.csv(tab, out_file, row.names = FALSE)
    cat("fan table:", out_file, "\n")
  }
  tab
}

## ================= RStudio: plot what is set at the top =================
# PLOT_VARIANT_NO_MAIN=1: only define plot_variant() (used by 09_results_figures.R, section 15)
NO_MAIN <- Sys.getenv("PLOT_VARIANT_NO_MAIN") == "1"
if (interactive() && !NO_MAIN) {
  print(plot_variant(GENE, CELLTYPE, VARIANT))
}

## ================= command line (Rscript) =================
if (!interactive() && !NO_MAIN) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) stop("usage: Rscript 10_plot_variant.R GENE CELLTYPE [VARIANT] | --file list.csv | --top N | --fan-table")
  if (args[1] == "--fan-table") {
    invisible(fan_table())
  } else if (args[1] %in% c("--file", "--top")) {
    todo <- if (args[1] == "--file") {
      read.csv(args[2])
    } else {
      hits %>% group_by(celltype) %>% arrange(p_gene) %>% slice_head(n = as.integer(args[2])) %>% ungroup() %>%
        select(gene_name, celltype) %>% as.data.frame()
    }
    if (!"variant_id" %in% names(todo)) todo$variant_id <- NA
    pdf_file <- file.path(plot_dir, paste0("variants_", if (args[1] == "--top") paste0("top", args[2]) else
      sub("\\.csv$", "", basename(args[2])), ".pdf"))
    pdf(pdf_file, width = 7, height = 6)
    for (r in seq_len(nrow(todo))) {
      p <- tryCatch(plot_variant(todo$gene_name[r], todo$celltype[r], todo$variant_id[r]),
                    error = function(e) { message("skipped ", todo$gene_name[r], ": ", conditionMessage(e)); NULL })
      if (!is.null(p)) print(p)
    }
    invisible(dev.off())
    cat("all pages:", pdf_file, "\n")
  } else {
    invisible(plot_variant(args[1], args[2], if (length(args) >= 3) args[3] else NA))
  }
}
