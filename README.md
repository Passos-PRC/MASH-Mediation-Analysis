# Bayesian Causal Mediation Analysis - MASH Drug Trials

This repository contains the analysis pipeline for a Bayesian causal mediation analysis of MASH/NASH clinical trials, investigating whether weight loss accounts for the histologic benefits of antifibrotic and metabolic therapies. It constructs treatment–placebo contrasts from extracted trial data and fits `brms`-based Bayesian random-effects models to decompose total effects into indirect (weight-mediated) and direct (weight-independent) components. Bayes Factors via the Savage–Dickey ratio are computed for each path.

---

## Repository structure

```
data/
  README.md                      column layout and SD harmonisation rules
  raw/                           place extraction spreadsheets here (git-ignored)

R/
  01_packages.R                  install and load all dependencies
  02_functions.R                 all analytic and plotting helpers
  03_example_weight_po1.R        end-to-end example: weight mediating PO1

output/                          generated files (git-ignored)
```

---

## How to run

1. **Install system dependency**: Stan requires a working C++ toolchain. See [mc-stan.org](https://mc-stan.org/users/interfaces/rstan).

2. **Place your data** in `data/raw/` following the column layout in `data/README.md`, or use the synthetic template.

3. **Set paths** via environment variables before running (or edit the defaults in the script):

```r
Sys.setenv(
  MEDIATOR_EXTRACTION_PATH = "data/raw/MEDIATORS_DATA_EXTRACTION.xlsx",
  MEDIATOR_OUTPUT_DIR      = "output"
)
```

4. **Run the example**:

```r
source("R/03_example_weight_po1.R")
```

This will produce `output/FOREST_WEIGHT_PO1.pdf`, `output/BUBBLE_WEIGHT_PO1.pdf`, and `output/mediation_weight_po1.csv`. The fitted Stan model is cached to disk so re-runs are instant.

---

## Model description

### Mediation framework

The analysis follows the product-of-coefficients approach:

- **a-path**: treatment / mediator effect, estimated per subgroup via a Normal–Normal conjugate Bayesian random-effects model using the within-study mediator contrast (active minus placebo).
- **b-path**: mediator / outcome slope, estimated globally via `brms` MCMC with a within–between decomposition (`wt_within` = arm deviation from study mean; `wt_mean` = study-level mediator mean) to avoid ecological confounding.
- **Indirect effect**: posterior product `a × b`, propagated via Monte Carlo draws.
- **Direct effect**: total effect − indirect effect.
- **Proportion mediated**: indirect / total, clipped to [0, 1] per draw.

### b-path formula

```
rd_resp | se(sqrt(var_rd_resp)) ~ wt_within + drug_class + wt_mean + (1 | study_id)
```

`rd_resp` is the arm-level absolute risk difference for the histologic endpoint; `var_rd_resp` is its within-arm variance (treated as known). The random intercept captures residual between-study heterogeneity.

### Priors

| Parameter | Prior |
|-----------|-------|
| Regression slopes | Normal(0, 0.5) |
| Intercept | Normal(0, 0.5) |
| Random-effect SD | Normal(0, 0.25) |

### Bayes Factors

BF₁₀ values are computed via the Savage–Dickey density ratio using a KDE of the posterior draws evaluated at zero. For indirect effects, a product-of-normals prior is used as the reference.

---

## Outputs

| File | Description |
|------|-------------|
| `FOREST_WEIGHT_PO1.pdf` | Forest plot + proportion-mediated panel |
| `BUBBLE_WEIGHT_PO1.pdf` | Trial-level scatter (weight loss vs. response) |
| `mediation_weight_po1.csv` | Full mediation table with 95% CrI and BF₁₀ |

---

## Citation

Please cite the final paper when using this code. If you adapt it for derivative analyses, consider linking back to this repository.

---
