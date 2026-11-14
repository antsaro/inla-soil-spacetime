# Space–time imputation of sparse soil records on an irrigated sugarcane estate

**Comparison of 25 latent Gaussian space–time structures (Besag, BYM2, AR(1), Knorr-Held interactions, Kronecker models) fitted with R-INLA, with hold-out validation, for 22 soil properties in Senegal**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R 4.1+](https://img.shields.io/badge/R-4.1%2B-276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![R-INLA](https://img.shields.io/badge/R--INLA-latent%20Gaussian%20models-6A4C93.svg)](https://www.r-inla.org/)
[![Quarto](https://img.shields.io/badge/Quarto-report-75AADB.svg?logo=quarto&logoColor=white)](https://quarto.org/)
[![sf](https://img.shields.io/badge/sf-spatial%20data-2E8B57.svg)](https://r-spatial.github.io/sf/)
[![ggplot2](https://img.shields.io/badge/ggplot2-figures-E36F47.svg)](https://ggplot2.tidyverse.org/)
[![Compositional data: ILR](https://img.shields.io/badge/texture-ILR%20coordinates-A0663A.svg)](#what-is-compared)
[![Study area: Senegal](https://img.shields.io/badge/study%20area-Senegal-00853F.svg)](#input-data)
[![Status: in development](https://img.shields.io/badge/status-in%20development-lightgrey.svg)](#running-it)

The data this was built for come from an irrigated sugarcane estate: about 850 fields followed over 23 seasons, with soil analysed only now and then (roughly once per replanting cycle). Most field × season cells are empty, a few fields were never sampled, and the analyses we want to run afterwards need a complete soil table. The question here is simply: *which space–time structure fills those gaps best?*

For every soil variable in the dataset, the code fits 25 structures with INLA and scores them on records that were left out of the fit.

## What is compared

All models share a Gaussian likelihood and fixed effects for management block and irrigation zone. They differ in the latent space–time part:

| Family | Structures |
|---|---|
| none | fixed effects only |
| space only | field iid, Besag (ICAR), BYM2 |
| time only | AR(1), RW1, RW2 on seasons |
| space + time | iid / Besag / BYM2 field effect + AR(1) seasons |
| space + time + interaction | Knorr-Held types I–IV |
| Kronecker space × time | Besag ⊗ AR(1), Besag ⊗ RW1, Besag ⊗ exchangeable, and variants |
| field-specific dynamics | AR(1) per field, linear trend per field |

The full list (S00–S24) is in section 3 of the Quarto document.

Neighbours are fields whose boundaries are less than 110 m apart. Priors are PC priors throughout, and responses are standardised before fitting.

**Responses.** Every numeric soil column is modelled. Strictly positive, right-skewed variables (skewness > 1) are log-transformed. Texture is never modelled as raw percentages: it goes in as pivot ILR coordinates, both for sand/silt/clay (`ILR3_1`, `ILR3_2`) and for the five fractions (`ILR5_1` to `ILR5_4`). The table of what was modelled, and how, is written to `responses_modelled.csv`.

**Scores.**
- In-sample: DIC, WAIC, mean log-CPO, marginal likelihood.
- Out-of-sample: RMSE, MAE, R², 95 % coverage, interval width, CRPS and log score, under two hold-out schemes:
  - `row`: 10 % of the measured field-seasons are removed at random. This is the gap-filling situation.
  - `field`: 10 % of the sampled fields are removed with all their records. This is the "field never sampled" situation.

Both schemes are repeated with different seeds (`N_REP`).

## Repository layout

```
analysis/model_comparison.qmd   main document: settings, models, parallel run, tables and figures
R/run_model_comparison.R        same pipeline as a plain R script (no report)
R/quick_screen.R                small sequential version for testing a few models on a few variables
R/install_packages.R            installs the packages, including INLA from its own repository
data/                           put the input .gpkg here (not tracked)
outputs/                        everything the code writes (not tracked)
```

## Requirements

- R ≥ 4.1
- INLA (from the INLA repository, not CRAN), sf, Matrix, dplyr, tidyr, ggplot2
- Quarto, if you want the rendered report rather than the plain script

```r
source("R/install_packages.R")
```

INLA is fast, but some of these models are not small. The Kronecker and Knorr-Held type III/IV structures build a latent field of fields × seasons (about 19,000 nodes) and need several GB of RAM per fit. Check memory before raising the number of workers.

## Input data

One GeoPackage, one row per field × season, with the soil columns empty when nothing was measured. The expected columns are described in [`data/README.md`](data/README.md). The data are not distributed with the code.

## Running it

**Quick check that everything works (a few minutes):**

```bash
VARS=pH MODELS="S02|S08|S14" SCHEMES=row Rscript R/quick_screen.R
```

**Full comparison.** Open `analysis/model_comparison.qmd`, check the settings block at the top, and render it:

```bash
quarto render analysis/model_comparison.qmd
```

or, without Quarto:

```bash
Rscript R/run_model_comparison.R
```

Points to know before a full run:

- **Parallel set-up.** The run is spread over `N_WORKERS` R processes with `parallel::parLapplyLB`. Each INLA call uses `THREADS_PER_FIT` threads. By default this uses (cores − 1) / 2 workers with 2 threads each.
- **It resumes.** Every finished variable × structure fit is saved in `outputs/model_comparison/fits/`. If the run stops, render again and only the missing fits are computed.
- **It is long.** With all variables, 25 structures and 3 repeats, expect a day or more on a workstation; the slow structures dominate. For a first pass, set `RUN_HEAVY <- FALSE` (the slow structures then run only on `HEAVY_VARS`) and/or `N_REP <- 1`.
- **Screening, then refit.** Screening uses `INT_STRATEGY <- "eb"` (hyperparameters at their mode). For the final numbers, refit the two or three best structures with `"ccd"` into a separate output folder.
- **INLA versions.** It was developed with an older INLA (19.09) and should run on current versions, but please open an issue if a model formula breaks.

## Outputs

In `outputs/model_comparison/`:

| File | Content |
|---|---|
| `all_results.csv` | one row per response × fixed-effect set × structure, with all criteria and scores |
| `ranking_overall.csv` | structures ranked by hold-out R², with the mean gain over the Besag ⊗ AR(1) model |
| `best_model_per_variable.csv` | best structure for each response |
| `holdout_predictions.csv` | every held-out record with its prediction and 95 % interval |
| `hyperparameters.csv` | posterior summaries of the hyperparameters |
| `responses_modelled.csv` | responses, transformation and number of measurements |
| `failed_fits.csv` | fits that did not converge or hit the time limit |
| `figures/` | heat map of hold-out R², gain over the reference model, row vs field validation, WAIC vs hold-out skill, observed vs predicted, calibration, effect of the fixed effects, best structure per response (PNG + PDF) |

## A few things worth knowing

- **WAIC and DIC are not hold-out skill.** Structures with an iid field × season term often get the best WAIC without predicting held-out records any better. That is why the ranking uses hold-out R², and why the comparison is shown explicitly in the WAIC figure.
- **The two hold-out schemes do not rank the models the same way.** Use the one that matches how the imputed values will be used.
- **Identifiability in types II–IV.** For Knorr-Held types II–IV only sum-to-zero constraints are set. The fitted values are fine, but the fixed effects of those models should not be interpreted.

## License

MIT, see [LICENSE](LICENSE).
