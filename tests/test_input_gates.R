source(file.path("R", "00_validate.R"))
source(file.path("R", "02_inference.R"))

# Relative abundance must be accepted as relative abundance.
ra <- matrix(c(0.2,0.8,0.4,0.6), nrow=2,
             dimnames=list(c("a","b"), c("s1","s2")))
x <- validate_abundance(ra, "relative_abundance")
stopifnot(x$data_type == "relative_abundance")

# Count-only wrappers must reject relative abundance before package invocation.
err <- try(require_genuine_counts(x), silent=TRUE)
stopifnot(inherits(err, "try-error"))

# Genuine integer counts pass the gate.
cnt <- matrix(c(2,8,4,6), nrow=2,
              dimnames=list(c("a","b"), c("s1","s2")))
y <- validate_abundance(cnt, "counts")
stopifnot(y$data_type == "counts")
require_genuine_counts(y)

# CLR coordinates are allowed to be negative and must reach inference.
source(file.path("R", "01_composition.R"))
source(file.path("R", "03_robustness.R"))
clr <- clr_transform(as_proportions(x), zero_method = "fixed", fixed = 1e-6)
stopifnot(any(clr < 0))
inf <- run_clr_inference(clr, c("vaginal", "c_section"),
                         c("vaginal", "c_section"), permutation_B = 9L)
stopifnot(length(inf) == 3L)

message("input gate tests passed")


# CLR coordinates may contain negative values and must pass transformed-matrix validation.
clr_demo <- matrix(c(-1, 1, -0.5, 0.5), nrow=2,
                   dimnames=list(c("a","b"), c("s1","s2")))
assert_matrix(clr_demo, "clr_demo", allow_negative=TRUE)

message("CLR transformed-value gate test passed")
