steps <- c(
  "01_load_bbs.R",
  "02_primary_robustness.R",
  "02b_age_sensitivity.R",
  "03_context_adjustment.R",
  "04_longitudinal.R",
  "04b_estimand_comparison.R",
  "05_make_figures.R"
)

for (s in steps) {
  message("\n== infant_birth_mode/", s, " ==")
  source(file.path("case_study", "infant_birth_mode", s),
         local = new.env(parent = globalenv()))
}

writeLines(capture.output(sessionInfo()),
           file.path("results", "sessionInfo.txt"))
message("\nCase study complete.")
