suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
})

source(file.path("R", "00_validate.R"))
source(file.path("R", "01_composition.R"))

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

tse <- curatedMetagenomicData(
  "2021-03-31.ShaoY_2019.relative_abundance",
  dryrun = FALSE,
  rownames = "short"
)[[1]]

# Primary estimand: early neonatal birth-mode contrast.
# Restrict first to a prespecified 0-7 day window, then keep one sample per
# infant closest to day 4. This avoids comparing infants sampled at materially
# different developmental ages merely because their earliest "newborn" sample
# occurred on a different day.
keep <- !is.na(tse$born_method) &
        tse$born_method %in% c("vaginal", "c_section") &
        !is.na(tse$age_category) &
        tse$age_category == "newborn" &
        !is.na(tse$subject_id) &
        !is.na(tse$infant_age) &
        tse$infant_age >= 0 &
        tse$infant_age <= 7

tse <- tse[, keep]

target_day <- 4
ord <- order(tse$subject_id, abs(tse$infant_age - target_day),
             tse$infant_age, na.last = TRUE)
tse <- tse[, ord]
tse <- tse[, !duplicated(tse$subject_id)]

ab <- as.matrix(assay(tse))
input <- validate_abundance(ab, data_type = "relative_abundance")
prop <- as_proportions(input)

meta <- as.data.frame(colData(tse))
rownames(meta) <- colnames(prop)
meta$born_method <- factor(meta$born_method, levels = c("vaginal", "c_section"))

saveRDS(
  list(input = input, proportions = prop, metadata = meta),
  file.path("results", "bbs_primary_input.rds")
)

writeLines(
  c(
    "Dataset: ShaoY_2019 / Baby Biome Study",
    "Input type: relative abundance",
    paste("Independent infants:", ncol(prop)),
    paste("Features before analysis filter:", nrow(prop)),
    paste("Primary age window: day 0-7; one sample per infant closest to day 4"),
    paste("Vaginal:", sum(meta$born_method == "vaginal")),
    paste("C-section:", sum(meta$born_method == "c_section"))
  ),
  file.path("results", "bbs_primary_provenance.txt")
)
