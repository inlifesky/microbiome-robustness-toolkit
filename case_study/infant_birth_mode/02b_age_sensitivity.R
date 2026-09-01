source(file.path("R", "00_validate.R"))
source(file.path("R", "01_composition.R"))
source(file.path("R", "02_inference.R"))

d <- readRDS(file.path("results", "bbs_primary_input.rds"))
meta <- d$metadata

prop <- filter_features(d$proportions, prevalence = 0.10)
clr <- clr_transform(prop, zero_method = "fixed", fixed = 1e-6)

g <- factor(meta$born_method, levels = c("vaginal", "c_section"))
age <- as.numeric(meta$infant_age)

# A. Age-adjusted CLR linear model.
rows_lm <- lapply(rownames(clr), function(f) {
  df <- data.frame(
    y = clr[f, ],
    born_method = g,
    infant_age = age
  )
  fit <- lm(y ~ born_method + infant_age, data = df)
  co <- summary(fit)$coefficients
  term <- "born_methodc_section"
  data.frame(
    feature = f,
    effect_age_adjusted = unname(co[term, "Estimate"]),
    p_age_adjusted = unname(co[term, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
})
age_lm <- do.call(rbind, rows_lm)
age_lm$q_age_adjusted <- p.adjust(age_lm$p_age_adjusted, "BH")

# B. Age-stratified label permutation.
# Birth-mode labels are shuffled only within observed sampling day, preserving
# the day-4/day-7 age composition under the null.
stratified_perm_p <- function(y, group, strata, B = 4999L) {
  group <- factor(group, levels = c("vaginal", "c_section"))
  obs <- mean(y[group == "c_section"]) - mean(y[group == "vaginal"])

  null <- replicate(B, {
    gp <- group
    for (s in unique(strata)) {
      ii <- which(strata == s)
      gp[ii] <- sample(gp[ii], replace = FALSE)
    }
    mean(y[gp == "c_section"]) - mean(y[gp == "vaginal"])
  })

  (1 + sum(abs(null) >= abs(obs))) / (B + 1)
}

unadjusted_effect <- rowMeans(clr[, g == "c_section", drop = FALSE]) -
                     rowMeans(clr[, g == "vaginal", drop = FALSE])

p_perm <- apply(
  clr, 1, stratified_perm_p,
  group = g, strata = age, B = 4999L
)

age_perm <- data.frame(
  feature = rownames(clr),
  effect_unadjusted = unadjusted_effect,
  p_age_stratified_permutation = p_perm,
  q_age_stratified_permutation = p.adjust(p_perm, "BH"),
  stringsAsFactors = FALSE
)

# C. Day-4-only sensitivity, because day 4 is the modal target day and
# comparison then has exact sampling-age matching.
keep4 <- which(age == 4)
day4 <- NULL
if (length(unique(g[keep4])) == 2L) {
  clr4 <- clr[, keep4, drop = FALSE]
  g4 <- droplevels(g[keep4])
  eff4 <- rowMeans(clr4[, g4 == "c_section", drop = FALSE]) -
          rowMeans(clr4[, g4 == "vaginal", drop = FALSE])
  p4 <- apply(clr4, 1, safe_wilcox, g = g4)
  day4 <- data.frame(
    feature = rownames(clr4),
    n_day4 = length(keep4),
    effect_day4_only = eff4,
    p_day4_only = p4,
    q_day4_only = p.adjust(p4, "BH"),
    stringsAsFactors = FALSE
  )
}

out <- merge(age_lm, age_perm, by = "feature", all = TRUE)
if (!is.null(day4)) out <- merge(out, day4, by = "feature", all = TRUE)

out$sign_age_adjusted_vs_unadjusted <-
  sign(out$effect_age_adjusted) == sign(out$effect_unadjusted)

if ("effect_day4_only" %in% names(out)) {
  out$sign_day4_vs_unadjusted <-
    sign(out$effect_day4_only) == sign(out$effect_unadjusted)
}

write.csv(
  out,
  file.path("results", "bbs_primary_age_sensitivity.csv"),
  row.names = FALSE
)

age_tab <- as.data.frame(table(
  infant_age = meta$infant_age,
  born_method = meta$born_method
))
write.csv(
  age_tab,
  file.path("results", "bbs_primary_age_by_birth_mode.csv"),
  row.names = FALSE
)
