# Data provenance

The repository itself contains no raw sequencing reads or participant-level case-study data.

The included infant birth-mode example retrieves public, uniformly processed Baby Biome Study metagenomic profiles at runtime through the Bioconductor package `curatedMetagenomicData`.

## Case-study source

- Original cohort: Baby Biome Study / Shao et al. 2019
- Public dataset identifier used for the primary analysis: `2021-03-31.ShaoY_2019.relative_abundance`
- Public dataset identifier used for context and longitudinal analyses: `ShaoY_2019.relative_abundance`
- Data type: species-level relative abundance
- Access layer: `curatedMetagenomicData`

The exact external-data contract, including the scripts that use each resource, is recorded in the repository-level `DATA_MANIFEST.tsv`.

## Derived files

Running the case study creates local RDS, CSV and figure outputs. These are generated artifacts and are not bundled in the source release. `sessionInfo.txt` records the executing R environment.

## Interpretation

Relative abundance does not provide absolute microbial abundance. The workflow does not attempt to reconstruct raw counts or absolute abundance from these profiles.
