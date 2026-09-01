# Compare two different estimands without forcing agreement:
# (1) repeated cross-sectional birth-mode contrasts within age windows;
# (2) predictions from subject-level mixed models.
#
# Agreement and disagreement are both retained in the diagnostic output.

source(file.path("R", "00_validate.R"))
source(file.path("R", "01_composition.R"))

windowed_path <- file.path("results", "bbs_windowed_group_effects.csv")
mixed_path <- file.path("results", "bbs_subject_mixed_models.csv")

if (file.exists(windowed_path) && file.exists(mixed_path)) {
  w <- read.csv(windowed_path, stringsAsFactors = FALSE)
  m <- read.csv(mixed_path, stringsAsFactors = FALSE)

  target <- data.frame(
    window = c("d0_7", "d8_30", "d91_plus"),
    target_day = c(4, 14, 210)
  )

  # Mixed model used standardized log1p(age). Reconstruct the scaling from
  # the same longitudinal metadata used by the case study.
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

  z_age <- log1p(meta$infant_age)
  mu <- mean(z_age, na.rm = TRUE)
  sig <- sd(z_age, na.rm = TRUE)

  target$z <- (log1p(target$target_day) - mu) / sig

  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(target))) {
    wi <- w[w$window == target$window[i], ]
    mm <- m
    mm$mixed_predicted_group_effect <-
      mm$group_effect_at_reference_age + mm$group_by_age * target$z[i]

    z <- merge(
      wi[, c("feature", "window", "effect")],
      mm[, c("feature", "mixed_predicted_group_effect",
             "converged", "singular")],
      by = "feature",
      all = FALSE
    )
    z$target_day <- target$target_day[i]
    z$sign_agrees <-
      sign(z$effect) == sign(z$mixed_predicted_group_effect)
    k <- k + 1L
    rows[[k]] <- z
  }

  out <- do.call(rbind, rows)
  names(out)[names(out) == "effect"] <- "windowed_effect"
  write.csv(
    out,
    file.path("results", "bbs_windowed_vs_mixed_estimand_comparison.csv"),
    row.names = FALSE
  )

  summary <- aggregate(
    sign_agrees ~ target_day,
    data = out,
    FUN = function(x) c(
      n = length(x),
      agree = sum(x),
      disagree = sum(!x),
      fraction_agree = mean(x)
    )
  )
  write.csv(
    summary,
    file.path("results", "bbs_windowed_vs_mixed_summary.csv"),
    row.names = FALSE
  )
}
