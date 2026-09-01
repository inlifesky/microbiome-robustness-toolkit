combine_inference_results <- function(result_list,
                                      q_threshold = 0.05) {
  methods <- names(result_list)
  if (!length(methods)) stop("result_list is empty.")

  features <- Reduce(union, lapply(result_list, function(x) x$feature))
  out <- data.frame(feature = features, stringsAsFactors = FALSE)

  for (m in methods) {
    x <- result_list[[m]]
    ii <- match(features, x$feature)
    out[[paste0(m, "_effect")]] <- x$effect[ii]
    out[[paste0(m, "_q")]] <- x$q[ii]
    out[[paste0(m, "_sig")]] <- !is.na(x$q[ii]) & x$q[ii] < q_threshold
  }

  sig_cols <- paste0(methods, "_sig")
  eff_cols <- paste0(methods, "_effect")

  out$n_significant <- rowSums(out[, sig_cols, drop = FALSE])

  out$sign_stable_among_significant <- vapply(seq_len(nrow(out)), function(i) {
    sig <- as.logical(out[i, sig_cols, drop = TRUE])
    eff <- as.numeric(out[i, eff_cols, drop = TRUE])
    s <- sign(eff[sig & is.finite(eff)])
    s <- s[s != 0]
    length(s) >= 2 && length(unique(s)) == 1
  }, logical(1))

  # Descriptive only: these tests are correlated if based on same CLR matrix.
  out
}

run_preprocessing_sensitivity <- function(input,
                                          group,
                                          group_levels,
                                          prevalence_grid = c(0.05, 0.10, 0.20),
                                          zero_methods = c("half_global_min", "fixed"),
                                          fixed_grid = c(1e-6, 1e-5),
                                          permutation_B = 999L) {
  prop0 <- as_proportions(input)
  all_runs <- list()
  id <- 0L

  for (prev in prevalence_grid) {
    prop <- filter_features(prop0, prevalence = prev)

    for (zm in zero_methods) {
      fvals <- if (zm == "fixed") fixed_grid else NA_real_

      for (fv in fvals) {
        id <- id + 1L
        clr <- clr_transform(
          prop,
          zero_method = zm,
          fixed = if (is.na(fv)) 1e-6 else fv
        )
        inf <- run_clr_inference(
          clr, group = group,
          group_levels = group_levels,
          permutation_B = permutation_B
        )
        comb <- combine_inference_results(inf)

        comb$spec_id <- id
        comb$prevalence <- prev
        comb$zero_method <- zm
        comb$fixed_value <- fv
        all_runs[[id]] <- comb
      }
    }
  }

  long <- do.call(rbind, all_runs)

  # Feature-level stability across preprocessing specifications.
  # A feature that disappears under stricter prevalence filtering is not
  # treated as fully stable: estimability_fraction records that explicitly.
  total_specs <- length(unique(long$spec_id))
  split_feature <- split(long, long$feature)

  # Per-specification mean effect (the three CLR tests share the same effect
  # estimator, so do not triple-count it for magnitude/rank summaries).
  spec_effect <- do.call(rbind, lapply(split_feature, function(x) {
    ecols <- grep("_effect$", names(x), value = TRUE)
    do.call(rbind, lapply(split(x, x$spec_id), function(z) {
      vals <- as.numeric(z[1, ecols, drop = TRUE])
      vals <- vals[is.finite(vals)]
      data.frame(
        feature = z$feature[1],
        spec_id = z$spec_id[1],
        effect = if (length(vals)) mean(vals) else NA_real_
      )
    }))
  }))

  # Within each preprocessing specification, rank by absolute effect size.
  spec_effect$abs_rank <- ave(
    -abs(spec_effect$effect), spec_effect$spec_id,
    FUN = function(z) rank(z, ties.method = "average", na.last = "keep")
  )

  summary <- do.call(rbind, lapply(split_feature, function(x) {
    q_cols <- grep("_q$", names(x), value = TRUE)
    qs <- unlist(x[, q_cols, drop = FALSE], use.names = FALSE)
    qs <- qs[is.finite(qs)]

    se <- spec_effect[spec_effect$feature == x$feature, ]
    effects <- se$effect[is.finite(se$effect)]
    signs <- sign(effects)
    signs <- signs[signs != 0]
    ranks <- se$abs_rank[is.finite(se$abs_rank)]

    data.frame(
      feature = x$feature[1],
      n_specifications_estimable = length(unique(x$spec_id)),
      estimability_fraction = length(unique(x$spec_id)) / total_specs,
      sign_stability = if (length(signs))
        max(table(signs)) / length(signs) else NA_real_,
      significant_fraction_across_tests =
        if (length(qs)) mean(qs < 0.05) else NA_real_,
      median_abs_effect = if (length(effects))
        median(abs(effects)) else NA_real_,
      median_abs_effect_rank = if (length(ranks))
        median(ranks) else NA_real_,
      effect_rank_iqr = if (length(ranks) > 1)
        IQR(ranks) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  list(all = long, specification_effects = spec_effect, summary = summary)
}
