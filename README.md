# HSI Directionality Analysis

Analysis code for the manuscript:

**Shifts in Clinician Management Decisions Following Exposure to an AI Model Validated to Assess Hydronephrosis Severity: A Simulated Pre-Deployment Study**

This repository contains the code used to prepare the analytic dataset, reproduce descriptive results, fit the mixed-effects models, generate the main dose-response figure, and perform exploratory subgroup and sensitivity analyses.

## Study analysis structure

The primary analysis evaluates how clinician management recommendations changed after exposure to the Hydronephrosis Severity Index (HSI).

Management intensity is ordered as:

1. Discharge
2. Repeat ultrasound
3. Diuretic renogram
4. Surgical referral

The primary directional outcomes are:

- **Escalation:** post-HSI recommendation is higher intensity than the initial recommendation
- **De-escalation:** post-HSI recommendation is lower intensity than the initial recommendation
- **No change:** post-HSI recommendation equals the initial recommendation

The HSI continuous predictor is analyzed as the **displayed HSI-predicted probability of future surgical intervention**, scaled per 10-percentage-point increase.

## Data

Patient-level and clinician-level source data are **not included in this repository**.

Place the source workbook locally at:

```text
data/Silent Trial Data.xlsx
```

or set the environment variable:

```bash
export HSI_DATA_PATH="/absolute/path/to/Silent Trial Data.xlsx"
```

Expected workbook sheets:

- `Decision Making Data`
- `Patient Variables`
- `Demographic`

The scripts contain explicit validation checks for the expected study structure:

- 23 clinicians
- 293 eligible case presentations per clinician after exclusion of case 267
- 6,739 eligible clinician-presentation observations
- 6,733 complete paired observations for management-transition analyses

If these checks fail, the scripts stop rather than silently proceeding.

## Repository layout

```text
analysis/
  00_setup.R
  01_prepare_data.R
  02_descriptive_tables.R
  03_mixed_models.R
  04_figure2.R
  05_subgroup_analysis.R
  06_exploratory_review_order.R
  07_ordinal_sensitivity.R

extras/
  optional_sankey.R

data/
  README.md

derived/
  # generated locally; not tracked

output/
  tables/
  figures/
  # generated locally; not tracked

run_all.R
```

## Required R packages

Core:

```r
install.packages(c(
  "readxl",
  "dplyr",
  "tidyr",
  "janitor",
  "lme4",
  "emmeans",
  "broom.mixed",
  "ggplot2",
  "patchwork",
  "openxlsx",
  "DescTools",
  "scales",
  "tibble",
  "here"
))
```

Sensitivity / optional:

```r
install.packages(c(
  "ordinal",
  "ggalluvial"
))
```

`MASS` and `splines` are included with standard R installations.

## Running the analysis

From the repository root:

```r
source("run_all.R")
```

This runs the core manuscript analyses:

1. data preparation
2. descriptive tables
3. primary mixed-effects models
4. Figure 2
5. subgroup analyses

Exploratory review-order and ordinal sensitivity analyses are intentionally not run by default. Run them separately if needed:

```r
source("analysis/06_exploratory_review_order.R")
source("analysis/07_ordinal_sensitivity.R")
```

## Important analysis notes

### Displayed green-risk probabilities

The study interface displayed green-risk HSI values as `0.0689` when the underlying model probability was below `0.069`. Analyses using the continuous HSI value therefore use the **displayed probability**, consistent with the clinician-facing interface.

### HSI risk thresholds

Prespecified thresholds are:

- Green: `< 0.069`
- Yellow: `>= 0.069 and < 0.34`
- Red: `>= 0.34`

### Eligible populations for directional models

To avoid impossible transitions:

- escalation: initial recommendation < surgical referral
- de-escalation: initial recommendation > discharge
- transition to surgical referral: initial recommendation is repeat ultrasound or diuretic renogram
- transition to discharge: initial recommendation is any non-discharge recommendation

For transition to discharge, all observed discharge transitions arose from an initial recommendation of repeat ultrasound. The continuous discharge analysis is therefore reported as a univariable mixed-effects model; adding initial management recommendation produced separation and unstable model estimation. Categorical adjusted-probability models are limited to the primary directional outcomes (escalation and de-escalation).

### Random effects

Mixed-effects models include random intercepts for:

- clinician
- case presentation

Adjusted models additionally include the initial management recommendation where applicable.

### Subgroup analyses

Specialty and experience analyses are exploratory.

Experience-level models additionally adjust for clinician specialty because specialty and training level are not distributed evenly across the sample. Fully stratified specialty-by-experience analyses are not used because some cells are sparse or empty.

### Review-order analysis

The source workbook used for the analysis does not encode exact session boundaries. The exploratory script therefore uses seven approximately equal **review-order blocks** as a proxy for the seven study sessions. This is explicitly labeled as exploratory and should not be described as exact session assignment unless the true session mapping is substituted.

## Reproducibility

The scripts do not automatically install packages and do not use `file.choose()`. Input and output locations are project-relative.

For a frozen package environment, initialize `renv` after confirming the analysis runs correctly:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

After reproducing the final manuscript outputs, save the R environment details with:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

## Code availability

If this repository is made public for publication, cite the permanent archived release (for example, a Zenodo DOI) in the manuscript rather than relying only on a mutable GitHub branch.
