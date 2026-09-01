assert_matrix <- function(x, name = "abundance", allow_negative = FALSE) {
  if (!is.matrix(x)) stop(name, " must be a matrix.")
  if (is.null(rownames(x)) || is.null(colnames(x)))
    stop(name, " must have feature row names and sample column names.")
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  if (!allow_negative && any(x < 0)) stop(name, " contains negative values.")
  invisible(TRUE)
}

validate_abundance <- function(x,
                               data_type = c("relative_abundance", "counts"),
                               tolerance = 0.05) {
  data_type <- match.arg(data_type)
  assert_matrix(x)

  if (data_type == "counts") {
    if (any(abs(x - round(x)) > sqrt(.Machine$double.eps)))
      stop("Count data must contain non-negative integer observations.")
    if (any(colSums(x) <= 0))
      stop("Every count sample must have a positive library size.")
  }

  if (data_type == "relative_abundance") {
    cs <- colSums(x)
    med <- median(cs[cs > 0])
    if (!is.finite(med)) stop("No positive sample totals.")
    # Permit proportions summing near 1 or percentages summing near 100.
    ok1 <- abs(cs - 1) <= tolerance
    ok100 <- abs(cs - 100) <= 100 * tolerance
    if (mean(ok1 | ok100) < 0.9)
      warning("Many sample totals are not close to 1 or 100; verify the input scale.")
  }

  structure(
    list(data = x, data_type = data_type),
    class = "microbiome_input"
  )
}

match_metadata <- function(input, metadata, sample_id_col = NULL) {
  stopifnot(inherits(input, "microbiome_input"))
  if (!is.data.frame(metadata)) stop("metadata must be a data.frame.")

  if (!is.null(sample_id_col)) {
    if (!sample_id_col %in% names(metadata)) stop("sample_id_col not found.")
    rownames(metadata) <- as.character(metadata[[sample_id_col]])
  }

  ids <- colnames(input$data)
  if (!all(ids %in% rownames(metadata)))
    stop("Some abundance sample IDs are missing from metadata.")

  metadata[ids, , drop = FALSE]
}

audit_subject_structure <- function(metadata, subject_col) {
  if (!subject_col %in% names(metadata)) stop("subject_col not found.")
  ids <- metadata[[subject_col]]
  tab <- table(ids, useNA = "ifany")
  list(
    n_samples = nrow(metadata),
    n_subjects = length(unique(ids[!is.na(ids)])),
    subjects_with_repeats = sum(tab > 1),
    max_samples_per_subject = if (length(tab)) max(tab) else NA_integer_
  )
}
