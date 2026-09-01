# Infant case-study diagnostics

This compact summary shows the diagnostics produced by the included public-data example.

## Primary sampling

- 554 infants: 296 vaginal birth and 258 caesarean birth.
- Samples in the prespecified day-0–7 analysis occurred on day 4 or day 7.
- Age-adjusted modelling produced no direction reversals relative to the reference analysis; the median absolute change in effect was approximately 0.027 CLR units.
- In the day-4-only sensitivity analysis, two features changed direction, but neither was FDR-significant.

## Preprocessing sensitivity

- 64 features entered the prevalence-sensitivity evaluation.
- Feature estimability varied across preprocessing specifications.
- Two features showed direction instability across the tested preprocessing grid.
- Zero-replacement strategy did not reverse feature direction in this case study, but it changed some CLR effect magnitudes enough to remain an important diagnostic.
- Effect-magnitude ranks were not perfectly fixed across specifications; rank stability is therefore reported separately from direction stability.

## Repeated observations

- Longitudinal input: 1,467 samples from 579 infants.
- 492 infants contributed repeated observations.
- All 61 fitted feature-level mixed models converged, were non-singular and had no substantive fitting warnings in this case study.
- Observed age windows contained 827 samples at day 0–7, 349 at day 8–30, none at day 31–90 and 291 at day 91+.

Windowed and mixed-model estimates were compared as separate estimands. Direction agreed for 57/61 features around day 4, 59/61 around day 14 and 50/61 around day 210.

## Context adjustment

For the 41 reported features in the neonatal context analysis, adding the available infant antibiotic-use context did not reverse the birth-mode coefficient direction. Absolute effect sizes decreased for some features and increased for others after adjustment; coefficient change is therefore reported directly.
