source(file.path("R", "00_validate.R"))
source(file.path("R", "01_composition.R"))
source(file.path("R", "05_context.R"))

# Pull full public metadata + profile for neonatal context analysis.
suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(dplyr)
})

meta_all <- sampleMetadata |>
  filter(study_name == "ShaoY_2019") |>
  filter(!is.na(infant_age), infant_age >= 0, infant_age <= 7) |>
  filter(born_method %in% c("vaginal", "c_section")) |>
  filter(antibiotics_current_use %in% c("no", "yes")) |>
  filter(!is.na(subject_id))

tse <- curatedMetagenomicData(
  "ShaoY_2019.relative_abundance",
  dryrun = FALSE,
  counts = FALSE,
  rownames = "short"
)[[1]]

ids <- intersect(meta_all$sample_id, colnames(tse))
meta <- meta_all[match(ids, meta_all$sample_id), , drop = FALSE]
ab <- as.matrix(assay(tse)[, ids, drop = FALSE])

# One neonatal sample per infant, closest to day 4.
ord <- order(meta$subject_id, abs(meta$infant_age - 4), meta$infant_age)
ids_ord <- seq_len(nrow(meta))[ord]
sel <- ids_ord[!duplicated(meta$subject_id[ord])]
meta <- meta[sel, , drop = FALSE]
ab <- ab[, sel, drop = FALSE]

input <- validate_abundance(ab, "relative_abundance")
prop <- filter_features(as_proportions(input), prevalence = 0.10)
clr <- clr_transform(prop, zero_method = "fixed", fixed = 1e-6)

rownames(meta) <- colnames(clr)
meta$born_method <- factor(meta$born_method, levels = c("vaginal", "c_section"))
meta$antibiotics_current_use <- factor(
  meta$antibiotics_current_use, levels = c("no", "yes")
)

ctx <- context_adjustment_sensitivity(
  clr = clr,
  metadata = meta,
  target_col = "born_method",
  target_levels = c("vaginal", "c_section"),
  context_cols = "antibiotics_current_use",
  extra_cols = "infant_age"
)

write.csv(ctx,
          file.path("results", "bbs_antibiotic_context_sensitivity.csv"),
          row.names = FALSE)
