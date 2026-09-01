suppressPackageStartupMessages({
  library(ggplot2)
})

sens <- read.csv(file.path("results", "bbs_preprocessing_sensitivity_summary.csv"))
ref <- read.csv(file.path("results", "bbs_reference_inference_stability.csv"))
lng <- read.csv(file.path("results", "bbs_windowed_group_effects.csv"))

# 1. Preprocessing sign stability
p1dat <- sens[order(-sens$estimability_fraction, -sens$sign_stability, -sens$median_abs_effect), ]
p1dat <- head(p1dat, 30)
p1dat$feature <- factor(p1dat$feature, levels = rev(p1dat$feature))
p1 <- ggplot(p1dat, aes(feature, sign_stability)) +
  geom_point(size = 2) +
  coord_flip() +
  ylim(0, 1) +
  labs(
    title = "Direction stability across preprocessing specifications",
    x = NULL, y = "Direction stability among estimable specifications"
  ) +
  theme_bw(base_size = 11)
ggsave(file.path("figures", "preprocessing_sign_stability.png"),
       p1, width = 7.5, height = 7, dpi = 220)

# 2. Inference stability (descriptive; correlated CLR tests)
p2dat <- ref[ref$n_significant > 0, ]
p2dat <- p2dat[order(p2dat$n_significant, p2dat$feature), ]
p2dat$feature <- factor(p2dat$feature, levels = p2dat$feature)
p2 <- ggplot(p2dat, aes(feature, n_significant,
                        shape = sign_stable_among_significant)) +
  geom_point(size = 2.5) +
  coord_flip() +
  scale_y_continuous(breaks = 1:3, limits = c(0.5, 3.2)) +
  labs(
    title = "Inference stability on a common CLR representation",
    subtitle = "Tests are complementary, not independent DA algorithms",
    x = NULL, y = "CLR tests with FDR < 0.05",
    shape = "Direction stable"
  ) +
  theme_bw(base_size = 11)
ggsave(file.path("figures", "clr_inference_stability.png"),
       p2, width = 7.5, height = max(4.5, 0.22*nrow(p2dat)+2), dpi = 220)

# 3. Windowed population contrasts.
# Do not connect d8-30 directly to d91+ because d31-90 has no observations.
stable <- head(
  sens$feature[
    order(-sens$estimability_fraction,
          -sens$sign_stability,
          -sens$median_abs_effect)
  ],
  12
)
p3dat <- lng[lng$feature %in% stable, ]
p3dat$window <- factor(
  p3dat$window,
  levels = c("d0_7", "d8_30", "d31_90", "d91_plus")
)

p3 <- ggplot(p3dat, aes(window, effect)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 1.8, alpha = 0.75) +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Birth-mode contrasts across observed infant-age windows",
    subtitle = "Points are not connected across the unobserved d31-90 window",
    x = "Age window",
    y = "CLR mean difference (C-section - vaginal)"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(
  file.path("figures", "windowed_birth_mode_contrast.png"),
  p3, width = 10, height = 8, dpi = 220
)
