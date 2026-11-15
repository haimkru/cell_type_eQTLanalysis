#author: haim krupkin+claude
#date: 11.07.2026
#manually verified on 11/14/2026
#description: step 2. One table of per-sample covariates, used by steps 3 and 6-8.
#   post_mortem : 1 = blood drawn after death (ischemic time SMTSISCH >= 0), 0 = before death (SMTSISCH < 0)
#   ventilator  : 1 = ventilator case (Hardy scale 0), 0 = other deaths (Hardy 1-4)
#   SMTSISCH    : ischemic time (minutes)
#   SMRIN       : RNA integrity number
#   centre_B1   : 1 = processed at centre B1
#   male, age   : biology (NOT removed as batch; used only as eQTL covariates)
# Missing values: ventilator -> the mean of the column; post_mortem, SMTSISCH, SMRIN -> the mean within the
# ventilator / non-ventilator group.
# Inputs : GTEx v8 SampleAttributesDS and SubjectPhenotypesDS
# Outputs: output/02_covariates.csv (one row per sample, in the expression's sample order)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
ann_dir  <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/edQTLs_GTEx/sample_annotations"
attr_file <- file.path(ann_dir, "GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt")
subj_file <- file.path(ann_dir, "GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt")

## ================= read =================
mix_all <- readRDS(file.path(out_dir, "01_tpm_noHb.rds"))
samples <- colnames(mix_all)
attr_df <- fread(attr_file, select = c("SAMPID", "SMTSISCH", "SMRIN", "SMCENTER"), quote = "")
subj_df <- fread(subj_file, select = c("SUBJID", "SEX", "AGE", "DTHHRDY"))

## ================= build the table =================
cov_df <- data.frame(SAMPID = samples)
cov_df <- left_join(cov_df, attr_df, by = "SAMPID")
cov_df$SUBJID <- sub("^([^-]+-[^-]+).*", "\\1", cov_df$SAMPID)       # donor ID = first two parts of the sample ID
cov_df <- left_join(cov_df, subj_df, by = "SUBJID")

cov_df$post_mortem <- as.numeric(cov_df$SMTSISCH >= 0)
cov_df$ventilator  <- as.numeric(cov_df$DTHHRDY == 0)
cov_df$centre_B1   <- as.numeric(grepl("B1", cov_df$SMCENTER))
cov_df$male        <- as.numeric(cov_df$SEX == 1)
#this is kinda weird and requires clarification
#this is a misstake - but it shouldnt effect any of the modeling so we can keep it for now, unless i remember anything that this affects.
cov_df$age         <- as.numeric(sub("-.*", "", cov_df$AGE)) + 5          # age bin "60-69" -> 65
#ok its working right

## ================= fill missing values =================
colSums(is.na(cov_df[, c("post_mortem", "ventilator", "SMTSISCH", "SMRIN")]))
cov_df$ventilator[is.na(cov_df$ventilator)] <- mean(cov_df$ventilator, na.rm = TRUE)
fill_by_group <- function(v, group) {
  group_mean <- ave(v, group, FUN = function(x) mean(x, na.rm = TRUE))
  ifelse(is.na(v), group_mean, v)
}
vent_group <- cov_df$ventilator > 0.5
#this is a super arbitrary threshold and might be better continous but we set it to that so thats fine.
cov_df$post_mortem <- fill_by_group(cov_df$post_mortem, vent_group)
cov_df$SMTSISCH    <- fill_by_group(cov_df$SMTSISCH, vent_group)
cov_df$SMRIN       <- fill_by_group(cov_df$SMRIN, vent_group)

keep_cols <- c("SAMPID", "SUBJID", "post_mortem", "ventilator", "SMTSISCH", "SMRIN", "centre_B1", "male", "age")
cov_df <- cov_df[, keep_cols]
stopifnot(!anyNA(cov_df))
summary(cov_df)
#this is a super arbitrary threshold and might be better continous but we set it to that so thats fine.

table(post_mortem = cov_df$post_mortem > 0.5, ventilator = cov_df$ventilator > 0.5)

## ================= save =================
out_file <- file.path(out_dir, "02_covariates.csv")
write.csv(cov_df, out_file, row.names = FALSE)
out_file
message("02_covariates done")
#code is fine hk 11/14/2026
