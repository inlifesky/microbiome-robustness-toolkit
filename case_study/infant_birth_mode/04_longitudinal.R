source(file.path("R", "00_validate.R"))
source(file.path("R", "01_composition.R"))
source(file.path("R", "02_inference.R"))
source(file.path("R", "04_longitudinal.R"))

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(dplyr)
})

meta_all <- sampleMetadata |>
  filter(study_name == "ShaoY_2019") |>
  filter(!is.na(infant_age), !is.na(subject_id)) |>
  filter(born_method %in% c("vaginal", "c_section"))

tse <- curatedMetagenomicData(
  "ShaoY_2019.relative_abundance",
  dryrun = FALSE,
  counts = FALSE,
  rownames = "short"
)[[1]]

ids <- intersect(meta_all$sample_id, colnames(tse))
meta <- meta_all[match(ids, meta_all$sample_id), , drop = FALSE]
ab <- as.matrix(assay(tse)[, ids, drop = FALSE])

input <- validate_abundance(ab, "relative_abundance")
prop <- filter_features(as_proportions(input), prevalence = 0.10)
clr <- clr_transform(prop, zero_method = "fixed", fixed = 1e-6)
rownames(meta) <- colnames(clr)

windowed <- windowed_group_effects(
  clr = clr,
  metadata = meta,
  subject_col = "subject_id",
  group_col = "born_method",
  age_col = "infant_age",
  group_levels = c("vaginal", "c_section"),
  breaks = c(-Inf, 7, 30, 90, Inf),
  labels = c("d0_7", "d8_30", "d31_90", "d91_plus"),
  target_ages = c(4, 14, 60, 210)
)

write.csv(windowed,
          file.path("results", "bbs_windowed_group_effects.csv"),
          row.names = FALSE)

# Mixed-model output is optional because lme4 may not be installed.
if (requireNamespace("lme4", quietly = TRUE)) {
  mixed <- fit_subject_mixed_models(
    clr = clr,
    metadata = meta,
    subject_col = "subject_id",
    group_col = "born_method",
    age_col = "infant_age",
    group_levels = c("vaginal", "c_section"),
    age_transform = "log1p",
    scale_age = TRUE
  )
  write.csv(mixed,
            file.path("results", "bbs_subject_mixed_models.csv"),
            row.names = FALSE)
}
