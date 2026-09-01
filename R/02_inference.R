safe_wilcox <- function(x, g) {
  tryCatch(stats::wilcox.test(x ~ g, exact = FALSE)$p.value,
           error = function(e) NA_real_)
}

safe_welch <- function(x, g) {
  tryCatch(stats::t.test(x ~ g, var.equal = FALSE)$p.value,
           error = function(e) NA_real_)
}

permutation_p <- function(x, g, B = 1999L) {
  g <- droplevels(factor(g))
  if (nlevels(g) != 2L) return(NA_real_)
  lev <- levels(g)
  obs <- mean(x[g == lev[2]]) - mean(x[g == lev[1]])

  null <- replicate(B, {
    gp <- sample(g, replace = FALSE)
    mean(x[gp == lev[2]]) - mean(x[gp == lev[1]])
  })

  (1 + sum(abs(null) >= abs(obs))) / (B + 1)
}

run_clr_inference <- function(clr,
                              group,
                              group_levels,
                              permutation_B = 1999L) {
  # CLR coordinates are centered log-ratios and therefore legitimately
  # contain both negative and positive values. The non-negative input gate
  # applies to abundance/count matrices, not transformed CLR coordinates.
  assert_matrix(clr, "clr", allow_negative = TRUE)
  if (length(group) != ncol(clr)) stop("group length does not match samples.")

  g <- factor(group, levels = group_levels)
  if (nlevels(droplevels(g)) != 2L) stop("Exactly two observed groups are required.")

  effect <- rowMeans(clr[, g == group_levels[2], drop = FALSE]) -
            rowMeans(clr[, g == group_levels[1], drop = FALSE])

  out <- list()

  p_w <- apply(clr, 1, safe_wilcox, g = g)
  out$clr_wilcoxon <- data.frame(
    feature = rownames(clr), effect = effect,
    p = p_w, q = p.adjust(p_w, "BH")
  )

  p_t <- apply(clr, 1, safe_welch, g = g)
  out$clr_welch <- data.frame(
    feature = rownames(clr), effect = effect,
    p = p_t, q = p.adjust(p_t, "BH")
  )

  p_perm <- apply(clr, 1, permutation_p, g = g, B = permutation_B)
  out$clr_permutation <- data.frame(
    feature = rownames(clr), effect = effect,
    p = p_perm, q = p.adjust(p_perm, "BH")
  )

  out
}

run_maaslin2_relative <- function(prop,
                                  metadata,
                                  fixed_effects,
                                  reference = NULL,
                                  output_dir,
                                  min_prevalence = 0) {
  if (!requireNamespace("Maaslin2", quietly = TRUE))
    stop("Maaslin2 is not installed.")

  if (!identical(colnames(prop), rownames(metadata)))
    stop("Metadata row order must match abundance columns.")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  Maaslin2::Maaslin2(
    input_data = as.data.frame(t(prop)),
    input_metadata = metadata,
    output = output_dir,
    fixed_effects = fixed_effects,
    reference = reference,
    normalization = "CLR",
    transform = "NONE",
    analysis_method = "LM",
    correction = "BH",
    min_prevalence = min_prevalence,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    cores = 1
  )
}

require_genuine_counts <- function(input) {
  if (!inherits(input, "microbiome_input")) stop("Use validate_abundance() first.")
  if (input$data_type != "counts")
    stop("This method requires genuine sequencing counts. Relative abundance cannot be converted to counts by scaling.")
  invisible(TRUE)
}

run_aldex2_counts <- function(input, group, ...) {
  require_genuine_counts(input)
  if (!requireNamespace("ALDEx2", quietly = TRUE)) stop("ALDEx2 is not installed.")
  ALDEx2::aldex(input$data, conditions = group, ...)
}

run_ancombc2_counts <- function(input, metadata, fix_formula, ...) {
  require_genuine_counts(input)
  if (!requireNamespace("ANCOMBC", quietly = TRUE)) stop("ANCOMBC is not installed.")

  ANCOMBC::ancombc2(
    data = input$data,
    taxa_are_rows = TRUE,
    meta_data = metadata,
    fix_formula = fix_formula,
    ...
  )
}
