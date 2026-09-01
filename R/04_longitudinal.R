pick_one_sample_per_subject_window <- function(metadata,
                                               subject_col,
                                               age_col,
                                               target_age) {
  subject <- metadata[[subject_col]]
  age <- metadata[[age_col]]
  ord <- order(subject, abs(age - target_age), age, na.last = TRUE)
  idx <- seq_len(nrow(metadata))[ord]
  idx[!duplicated(subject[ord])]
}

windowed_group_effects <- function(clr,
                                   metadata,
                                   subject_col,
                                   group_col,
                                   age_col,
                                   group_levels,
                                   breaks,
                                   labels,
                                   target_ages) {
  if (length(labels) != length(target_ages))
    stop("labels and target_ages must have the same length.")
  if (!identical(colnames(clr), rownames(metadata)))
    stop("metadata rows must match CLR sample columns.")

  metadata$.window <- cut(
    metadata[[age_col]],
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE
  )

  rows <- list()
  k <- 0L

  for (j in seq_along(labels)) {
    w <- labels[j]
    idx <- which(metadata$.window == w &
                 metadata[[group_col]] %in% group_levels &
                 !is.na(metadata[[subject_col]]) &
                 !is.na(metadata[[age_col]]))
    if (!length(idx)) next

    mw <- metadata[idx, , drop = FALSE]
    local <- pick_one_sample_per_subject_window(
      mw, subject_col, age_col, target_ages[j]
    )
    idx <- idx[local]
    mw <- metadata[idx, , drop = FALSE]
    cw <- clr[, idx, drop = FALSE]

    g <- factor(mw[[group_col]], levels = group_levels)

    for (f in rownames(cw)) {
      y <- cw[f, ]
      k <- k + 1L
      rows[[k]] <- data.frame(
        feature = f,
        window = as.character(w),
        n_group_a = sum(g == group_levels[1]),
        n_group_b = sum(g == group_levels[2]),
        effect = mean(y[g == group_levels[2]]) - mean(y[g == group_levels[1]]),
        p = safe_wilcox(y, g),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  out$q_within_window <- ave(out$p, out$window, FUN = function(z) p.adjust(z, "BH"))
  out
}

fit_subject_mixed_models <- function(clr,
                                     metadata,
                                     subject_col,
                                     group_col,
                                     age_col,
                                     group_levels,
                                     age_transform = c("log1p", "linear"),
                                     scale_age = TRUE,
                                     min_subjects_with_repeats = 20L) {
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("lme4 is required for mixed models.")
  if (!identical(colnames(clr), rownames(metadata)))
    stop("metadata rows must match CLR sample columns.")

  age_transform <- match.arg(age_transform)
  meta <- metadata
  meta$.group <- factor(meta[[group_col]], levels = group_levels)
  meta$.subject <- factor(meta[[subject_col]])

  raw_age <- meta[[age_col]]
  if (age_transform == "log1p") {
    if (any(raw_age < 0, na.rm = TRUE))
      stop("log1p age transform requires non-negative ages.")
    meta$.age <- log1p(raw_age)
  } else {
    meta$.age <- raw_age
  }

  if (scale_age) meta$.age <- as.numeric(scale(meta$.age))

  keep <- complete.cases(meta[, c(".group", ".subject", ".age")])
  meta <- meta[keep, , drop = FALSE]
  cc <- clr[, keep, drop = FALSE]

  repeat_tab <- table(meta$.subject)
  n_repeat_subjects <- sum(repeat_tab >= 2)
  if (n_repeat_subjects < min_subjects_with_repeats) {
    warning("Only ", n_repeat_subjects,
            " subjects have repeated observations; mixed-model estimates may be unstable.")
  }

  rows <- lapply(rownames(cc), function(f) {
    df <- meta
    df$.y <- cc[f, ]

    warns <- character()
    fit <- withCallingHandlers(
      tryCatch(
        lme4::lmer(.y ~ .group * .age + (1 | .subject),
                   data = df, REML = FALSE),
        error = function(e) NULL
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (is.null(fit)) return(NULL)

    co <- summary(fit)$coefficients
    interaction_row <- grep("^\\.group.*:\\.age$|^\\.age:\\.group",
                            rownames(co), value = TRUE)
    group_row <- grep("^\\.group", rownames(co), value = TRUE)
    group_row <- setdiff(group_row, interaction_row)

    opt_messages <- fit@optinfo$conv$lme4$messages
    converged <- is.null(opt_messages)
    singular <- lme4::isSingular(fit, tol = 1e-4)

    data.frame(
      feature = f,
      age_transform = age_transform,
      n_samples = nrow(df),
      n_subjects = length(unique(df$.subject)),
      n_subjects_with_repeats = n_repeat_subjects,
      group_effect_at_reference_age =
        if (length(group_row)) co[group_row[1], "Estimate"] else NA_real_,
      group_by_age =
        if (length(interaction_row)) co[interaction_row[1], "Estimate"] else NA_real_,
      group_by_age_se =
        if (length(interaction_row)) co[interaction_row[1], "Std. Error"] else NA_real_,
      converged = converged,
      singular = singular,
      warning_message = paste(unique(c(warns, opt_messages)), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
