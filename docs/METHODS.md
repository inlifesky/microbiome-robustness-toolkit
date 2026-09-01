# Methods

## 1. Start with the data type

The workflow distinguishes `counts` from `relative_abundance` before any differential-abundance method is called.

Genuine sequencing counts retain library-size information and can support count-based methods when their assumptions are otherwise appropriate. Relative-abundance tables are treated as compositions. Count-based methods require genuine count input.

## 2. Define the estimand before testing

Examples include:

- a cross-sectional group contrast in a prespecified age window;
- a covariate-adjusted group association;
- a group-by-age term in repeated observations.

The longitudinal module therefore separates population contrasts within age windows from subject-level mixed models.

## 3. Compositional representation

Relative-abundance tables are closed to proportions and transformed with centered log-ratios (CLR). Zero replacement is an explicit parameter because log-ratio transformation requires positive values.

Available strategies include a small fixed replacement, half of the global minimum positive abundance, and optional multiplicative replacement through `zCompositions`.

## 4. Inference stability

For a common CLR representation, the toolkit can compare:

- Wilcoxon rank-sum test;
- Welch two-sample test;
- label-permutation test.

These tests share the same CLR representation and are used to assess sensitivity to the inferential statistic.

`Maaslin2` can be used as a separate multivariable modelling layer. ALDEx2 and ANCOM-BC2 wrappers are available only for validated count input.

## 5. Preprocessing robustness

A sensitivity grid can vary:

- prevalence threshold;
- zero-replacement strategy;
- fixed replacement value.

For each feature the toolkit records:

- estimability across specifications;
- effect direction stability;
- fraction of tests meeting the FDR threshold;
- median absolute effect;
- effect-magnitude rank and rank IQR.

The outputs are kept as separate stability dimensions. A feature may have a stable direction while its rank changes, or it may be estimable only under permissive filtering.

## 6. Sampling-age sensitivity

In rapidly developing systems, selecting each subject's earliest available sample is not sufficient control for age.

The infant case study prespecifies day 0–7 and selects one sample per infant closest to day 4. It then checks the result with:

- an age-adjusted CLR linear model;
- permutation stratified by observed sampling day;
- a day-4-only sensitivity subset.

## 7. Longitudinal analyses

### Repeated cross-sectional windows

At most one sample per subject is retained in each window. The output is a group contrast for that age window.

### Subject-level mixed model

Repeated samples are represented with a subject random intercept. Age transformation is explicit, and outputs include convergence and singular-fit diagnostics.

The two approaches estimate different quantities, so their agreement and disagreement are reported directly.

Missing age windows are shown as missing. The plotting code does not draw a continuous line through an unobserved interval.

## 8. Context adjustment

The context module compares a target coefficient before and after adding prespecified covariates. It reports signed change, absolute-effect change and direction changes.

Changes in the adjusted coefficient are reported as model sensitivity to covariate adjustment.

## 9. Extending the toolkit

For a new dataset:

1. validate the input type and metadata alignment;
2. define the biological estimand and independent sampling unit;
3. choose preprocessing sensitivity ranges appropriate to the dataset;
4. run inference and effect-stability checks;
5. add covariates or repeated-measure models only when supported by the design;
6. carry forward signals whose stability profile is appropriate for the next biological question.
