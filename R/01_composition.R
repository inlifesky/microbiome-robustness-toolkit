as_proportions <- function(input) {
  stopifnot(inherits(input, "microbiome_input"))
  x <- input$data
  if (input$data_type == "counts") {
    x <- sweep(x, 2, colSums(x), "/")
  } else {
    cs <- colSums(x)
    if (median(cs[cs > 0]) > 2) x <- x / 100
    # Re-close to 1 to remove small numerical scale differences.
    x <- sweep(x, 2, colSums(x), "/")
  }
  x
}

filter_features <- function(prop,
                            prevalence = 0.10,
                            min_max_abundance = 1e-4) {
  assert_matrix(prop, "prop")
  keep <- rowMeans(prop > 0) >= prevalence &
          apply(prop, 1, max) >= min_max_abundance
  prop[keep, , drop = FALSE]
}

replace_zeros <- function(prop,
                          method = c("half_global_min", "fixed", "multiplicative"),
                          fixed = 1e-6) {
  method <- match.arg(method)
  assert_matrix(prop, "prop")

  if (!any(prop == 0)) return(prop)

  if (method == "half_global_min") {
    pos <- prop[prop > 0]
    if (!length(pos)) stop("No positive values available for replacement.")
    pc <- min(pos) / 2
    y <- prop
    y[y == 0] <- pc
    y <- sweep(y, 2, colSums(y), "/")
    return(y)
  }

  if (method == "fixed") {
    if (!is.numeric(fixed) || length(fixed) != 1 || fixed <= 0)
      stop("fixed must be a positive scalar.")
    y <- prop
    y[y == 0] <- fixed
    y <- sweep(y, 2, colSums(y), "/")
    return(y)
  }

  if (!requireNamespace("zCompositions", quietly = TRUE))
    stop("zCompositions is required for multiplicative replacement.")

  # cmultRepl expects samples in rows.
  y <- zCompositions::cmultRepl(t(prop), label = 0, method = "CZM")
  t(y)
}

clr_transform <- function(prop,
                          zero_method = "multiplicative",
                          fixed = 1e-6) {
  y <- replace_zeros(prop, method = zero_method, fixed = fixed)
  logy <- log(y)
  sweep(logy, 2, colMeans(logy), "-")
}
