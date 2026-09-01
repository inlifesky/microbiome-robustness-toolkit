context_adjustment_sensitivity <- function(clr,
                                           metadata,
                                           target_col,
                                           target_levels,
                                           context_cols,
                                           extra_cols = NULL) {
  if (!identical(colnames(clr), rownames(metadata)))
    stop("metadata rows must match CLR sample columns.")

  meta <- metadata
  meta$.target <- factor(meta[[target_col]], levels = target_levels)

  all_cov <- unique(c(context_cols, extra_cols))
  needed <- c(".target", all_cov)
  keep <- complete.cases(meta[, needed, drop = FALSE])

  meta <- meta[keep, , drop = FALSE]
  cc <- clr[, keep, drop = FALSE]

  f0 <- stats::as.formula(".y ~ .target")
  rhs <- paste(c(".target", all_cov), collapse = " + ")
  f1 <- stats::as.formula(paste(".y ~", rhs))

  rows <- lapply(rownames(cc), function(f) {
    df <- meta
    df$.y <- cc[f, ]

    m0 <- lm(f0, data = df)
    m1 <- lm(f1, data = df)

    term <- grep("^\\.target", names(coef(m0)), value = TRUE)[1]
    if (is.na(term)) return(NULL)

    b0 <- unname(coef(m0)[term])
    b1 <- unname(coef(m1)[term])

    p0 <- summary(m0)$coefficients[term, "Pr(>|t|)"]
    p1 <- if (term %in% rownames(summary(m1)$coefficients))
      summary(m1)$coefficients[term, "Pr(>|t|)"] else NA_real_

    data.frame(
      feature = f,
      unadjusted_effect = b0,
      adjusted_effect = b1,
      absolute_change = b1 - b0,
      proportional_abs_effect_change =
        if (is.finite(b0) && abs(b0) > 0)
          abs(b1) / abs(b0) - 1 else NA_real_,
      sign_changed = is.finite(b0) && is.finite(b1) && sign(b0) != sign(b1),
      unadjusted_p = p0,
      adjusted_p = p1,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$unadjusted_q <- p.adjust(out$unadjusted_p, "BH")
  out$adjusted_q <- p.adjust(out$adjusted_p, "BH")
  out
}
