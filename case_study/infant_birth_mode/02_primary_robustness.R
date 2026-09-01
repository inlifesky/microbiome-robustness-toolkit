source(file.path("R", "00_validate.R"))
source(file.path("R", "01_composition.R"))
source(file.path("R", "02_inference.R"))
source(file.path("R", "03_robustness.R"))

d <- readRDS(file.path("results", "bbs_primary_input.rds"))
input <- d$input
meta <- d$metadata

sens <- run_preprocessing_sensitivity(
  input = input,
  group = meta$born_method,
  group_levels = c("vaginal", "c_section"),
  prevalence_grid = c(0.05, 0.10, 0.20),
  zero_methods = c("half_global_min", "fixed"),
  fixed_grid = c(1e-6, 1e-5),
  permutation_B = 999L
)

write.csv(sens$all,
          file.path("results", "bbs_preprocessing_sensitivity_all.csv"),
          row.names = FALSE)
write.csv(sens$summary,
          file.path("results", "bbs_preprocessing_sensitivity_summary.csv"),
          row.names = FALSE)
write.csv(sens$specification_effects,
          file.path("results", "bbs_preprocessing_specification_effects.csv"),
          row.names = FALSE)

# One prespecified reference specification for readable reporting.
prop <- filter_features(d$proportions, prevalence = 0.10)
clr <- clr_transform(prop, zero_method = "fixed", fixed = 1e-6)

inf <- run_clr_inference(
  clr,
  group = meta$born_method,
  group_levels = c("vaginal", "c_section"),
  permutation_B = 1999L
)
combined <- combine_inference_results(inf)
write.csv(combined,
          file.path("results", "bbs_reference_inference_stability.csv"),
          row.names = FALSE)

saveRDS(list(prop = prop, clr = clr, inference = inf),
        file.path("results", "bbs_reference_analysis.rds"))
