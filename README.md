# Microbiome robustness toolkit

An R workflow for evaluating how microbiome associations change across analytical specifications.

Microbiome results can vary with data representation, zero handling, feature filtering, covariate adjustment, sampling age and repeated-measure modelling. The toolkit reports these sources of variation separately.

**[Project overview / GitHub Pages](https://inlifesky.github.io/microbiome-robustness-toolkit/)**

## What it does

The workflow is organised around six checks:

1. **Input integrity** — distinguishes genuine sequencing counts from relative-abundance tables and checks sample/metadata alignment.
2. **Compositional preprocessing** — CLR transformation with explicit zero-replacement choices.
3. **Inference stability** — complementary tests on the same CLR representation, plus optional multivariable modelling.
4. **Preprocessing stability** — prevalence, zero handling, effect direction, effect magnitude and rank stability.
5. **Longitudinal structure** — separates repeated cross-sectional contrasts from subject-level mixed models.
6. **Context adjustment** — compares the target association before and after covariate adjustment without treating coefficient change as causal mediation.

**Within-representation inference stability** and **cross-algorithm evidence** are reported separately. Count-based methods such as ALDEx2 and ANCOM-BC2 are available when genuine count data are supplied. Relative-abundance tables remain on the compositional analysis path.

## Typical use

A practical place for this workflow is between data preparation and biological interpretation:

```text
microbiome feature table
        ↓
input / metadata audit
        ↓
compositional + preprocessing sensitivity
        ↓
inference / covariate / longitudinal checks
        ↓
stable, unstable, or analysis-dependent signals
        ↓
biological modelling, validation, or follow-up
```

The same workflow can be used to re-analyse a published dataset or to screen a new dataset before deeper biological modelling.

## Included case study

`case_study/infant_birth_mode/` demonstrates the workflow with public Baby Biome Study shotgun-metagenomic relative-abundance profiles accessed through `curatedMetagenomicData`.

The example demonstrates the diagnostic workflow. The primary comparison is restricted to day 0–7, with one sample per infant selected closest to day 4. Age sensitivity, preprocessing sensitivity, context adjustment, windowed contrasts and mixed-model diagnostics are reported separately.

In the case study, the primary analysis includes 554 infants and the longitudinal analysis includes 1,467 samples from 579 infants. Effect direction is relatively stable for many features, while effect magnitude, estimability and rank vary across analytical specifications.

See [`example_outputs/CASE_STUDY_DIAGNOSTICS.md`](example_outputs/CASE_STUDY_DIAGNOSTICS.md) for the compact execution summary.

## Run the toolkit

Load the reusable modules:

```r
source("R/00_validate.R")
source("R/01_composition.R")
source("R/02_inference.R")
source("R/03_robustness.R")
source("R/04_longitudinal.R")
source("R/05_context.R")
```

Run the included case study from the repository root:

```bash
Rscript case_study/infant_birth_mode/run_case_study.R
```

The case study retrieves the public data at runtime. Dataset identifiers and data flow are recorded in [`DATA_MANIFEST.tsv`](DATA_MANIFEST.tsv). Generated outputs are written to `results/` and `figures/`, and the R environment is captured in `results/sessionInfo.txt`.

## Repository structure

```text
R/                         reusable analysis modules
case_study/                worked public-data example
tests/                     input-gate regression tests
docs/                      methods and data provenance
example_outputs/           compact case-study diagnostics
index.html                 lightweight project overview
DATA_MANIFEST.tsv          external-data contract
REFERENCES.md              data and software references
ACKNOWLEDGEMENTS.md         data/software acknowledgements
THIRD_PARTY_NOTICES.md      third-party software and data license notices
CITATION.cff               citation metadata
```

## Documentation

- [Methods](docs/METHODS.md)
- [Data provenance](docs/DATA_PROVENANCE.md)
- [References](REFERENCES.md)
- [Acknowledgements](ACKNOWLEDGEMENTS.md)
- [Third-party notices](THIRD_PARTY_NOTICES.md)
- [Project overview](index.html)

## Dependencies

Core R functions are supplemented, where required, by `dplyr`, `tidyr`, `ggplot2`, `Maaslin2`, `lme4`, `zCompositions`, `ALDEx2`, `ANCOMBC`, `curatedMetagenomicData` and `SummarizedExperiment`. Not every package is required for every analysis path.

The included case study has been run successfully on real public data. Exact package versions should be recorded from the generated `sessionInfo.txt` when reproducing or adapting the workflow.

## License

MIT. Data remain subject to the terms of their original repositories and publications.
