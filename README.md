# bcnma

**Bayesian Component Network Meta-Analysis with LASSO Interactions**

`bcnma` fits Bayesian random-effects component network meta-analysis (cNMA)
for binary outcomes using JAGS. It supports three interaction models —
additive, SSVS-LASSO (equiprobable), and SSVS-LASSO (informative) — and
provides `netmeta`-style output for users already familiar with the
`meta`/`netmeta` ecosystem.

The package was developed by packaging and generalising analysis scripts
originally written for the **CBT-I Component NMA** study:

> Furukawa Y, Sakata M, Yamamoto R, et al. (2024). Components and delivery
> formats of cognitive behavioral therapy for chronic insomnia in adults.
> *JAMA Psychiatry*, 81(2), 130–140.
> https://doi.org/10.1001/jamapsychiatry.2023.5060

The bundled example dataset (`inst/extdata/data_CBTICNMA.csv`) is the
publicly available dataset from that study
([github.com/ykfrkw/CBTICNMA](https://github.com/ykfrkw/CBTICNMA)).

---

## Installation

### Prerequisites

1. **JAGS** (Just Another Gibbs Sampler) must be installed separately:
   https://sourceforge.net/projects/mcmc-jags/

2. Install R package dependencies:

```r
install.packages(c("rjags", "MCMCvis", "ggplot2", "ggrepel"))
```

### Install bcnma from GitHub

```r
# Install remotes if needed
install.packages("remotes")

# Install bcnma
remotes::install_github("ykfrkw/bcnma")
```

---

## Quick start

```r
library(bcnma)
library(dplyr)

# Load the bundled CBTICNMA dataset (248 studies, 17 CBT-I components)
dat_raw <- read.csv(system.file("extdata", "data_CBTICNMA.csv", package = "bcnma"))
names(dat_raw)[1] <- "study"  # fix BOM character in first column

dat <- dat_raw |>
  select(study, treatment_component, r, n) |>
  group_by(study, treatment_component) |>
  summarise(r = sum(r), n = sum(n), .groups = "drop")

components <- c("se", "sd", "cr", "th", "cw", "sr", "sc", "re", "pi",
                "w",  "ns", "he", "tg", "ind", "gp", "ff", "ae")

# Fit LASSO model (all pairwise interactions with SSVS)
fit <- bcnma(
  data         = dat,
  studlab      = study,
  treat        = treatment_component,
  event        = r,
  n            = n,
  components   = components,
  interactions = "lasso",
  seed         = 42
)

print(fit)
forest(fit)
plot_interactions(fit)
```

---

## Input data format

Long format, **one row per arm per study**:

| study      | treatment_component | r  | n  |
|------------|--------------------|----|-----|
| Smith2001  | se+sd+cr           | 20 | 50 |
| Smith2001  | w                  | 10 | 50 |
| Jones2010  | se+sd+cr+th+sr+sc  | 30 | 60 |
| Jones2010  | se                 | 15 | 55 |

- **`treat` column**: component names joined by `component_separator` (default `"+"`)
- **`components` vector**: all component abbreviations; freely customisable for any number of components
- **`event` / `n`**: binary outcome numerator and denominator per arm

---

## Three interaction models

### `interactions = "none"` — Additive (no interactions)

Pure additive model. The effect of a combination equals the sum of its
component effects. Use as a baseline comparison.

```r
bcnma(..., interactions = "none")
```

### `interactions = "lasso"` — Equiprobable SSVS (recommended default)

All pairwise interactions are evaluated with equal prior inclusion
probability (`lasso_prob = 0.5`). The data determine which interactions
are supported.

```r
bcnma(..., interactions = "lasso", lasso_prob = 0.5, g = 100)
```

### `interactions = "lasso_informative"` — Informative SSVS

Only interactions among `active_components` are modelled stochastically;
all other interactions are fixed to zero. Use when prior knowledge or
clinical theory supports focusing on specific component pairs.

```r
bcnma(...,
  interactions      = "lasso_informative",
  active_components = c("se", "sd", "cr", "th", "sr", "sc"),
  lasso_prob        = 0.8)
```

---

## Key arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `components` | — | **Required.** Character vector of component abbreviations. No hard-coded limit. |
| `component_separator` | `"+"` | Separator in the `treat` column |
| `interactions` | `"lasso"` | `"none"`, `"lasso"`, or `"lasso_informative"` |
| `active_components` | `NULL` | Required for `"lasso_informative"` |
| `lasso_prob` | `0.5` | Prior P(interaction included) |
| `g` | `100` | Spike-slab variance ratio (slab = g × spike variance) |
| `tau_prior` | `c(-1.67, 1/1.472^2)` | Log-normal prior on τ (Turner et al. 2015) |
| `n.chains` | `3` | MCMC chains (4 recommended for LASSO models) |
| `n.adapt` | `2000` | JAGS adaptation iterations |
| `n.burnin` | `5000` | Burn-in iterations after adaptation |
| `n.iter` | `20000` | Posterior sampling iterations per chain |
| `seed` | `NULL` | Random seed for reproducibility |

---

## Output functions

### `print(fit)` / `summary(fit)`

Component effects table: log-OR and OR with 95% credible interval (CrI)
for both additive and LASSO models. Also reports the number of interactions
with posterior inclusion probability ≥ 0.5.

```r
print(fit)
print(fit, model = "lasso")     # LASSO model only
summary(fit)                    # full output with MCMC settings
```

### `forest(fit)`

Forest plot of component OR estimates. Displays both additive and LASSO
estimates side-by-side when `model = "both"` (default).

```r
forest(fit)
forest(fit, model = "lasso")   # LASSO model only
```

Returns a `ggplot` object invisibly for further customisation.

### `plot_interactions(fit)`

Scatter plot of median interaction effect (γ, x-axis, log-OR scale) vs.
posterior inclusion probability (y-axis) for all pairwise interactions.
Top `top_n` pairs by |γ| are labelled.

```r
plot_interactions(fit, top_n = 10, threshold = 0.5)
```

### `summarise_interactions(fit)`

Data frame of all interaction pairs sorted by |γ| (descending):

| Column | Description |
|--------|-------------|
| `component_pair` | Component names (e.g. `"cr - th"`) |
| `median_gamma` | Posterior median of γ (log-OR) |
| `OR_interaction` | exp(γ̂) |
| `inclusion_prob` | Posterior P(Ind = 1) |
| `selected` | `TRUE` if `inclusion_prob >= threshold` |

```r
tbl <- summarise_interactions(fit, threshold = 0.5)
print(tbl)
```

### `compare_combinations(fit, combo1, combo2)`

Posterior OR (and 95% CrI) for combination A vs. combination B, accounting
for all relevant interaction terms automatically.

```r
# Full CBT-I package vs. sleep hygiene only
compare_combinations(fit,
  combo1 = c("se", "sd", "cr", "th", "sr", "sc"),
  combo2 = c("se"))

# Any combination vs. no-component baseline
compare_combinations(fit,
  combo1 = c("cr", "th", "sr"),
  combo2 = character(0))

# Additive model comparison (no interaction adjustment)
compare_combinations(fit,
  combo1 = c("se", "sd", "cr", "th"),
  combo2 = c("se"),
  model  = "additive")
```

---

## How to read the results

### Component effects (from `print(fit)`)

Each row shows the **main effect** of one component vs. "no component
present" (baseline), pooled across all studies in the network.

| Column | Interpretation |
|--------|----------------|
| `log-OR [95% CrI]` | Log odds ratio; positive = beneficial |
| `OR [95% CrI]` | Odds ratio on the natural scale; OR > 1 = beneficial |
| `tau` | Between-study SD of log-ORs; measures statistical heterogeneity |

**Example from CBTICNMA data (additive model):**

```
  Comp      log-OR [95% CrI]          OR [95% CrI]
  cr        0.53 [ 0.21;  0.86]       1.71 [1.23; 2.36]  ← beneficial
  th        0.59 [ 0.25;  0.94]       1.81 [1.28; 2.57]  ← beneficial
  sr        0.47 [ 0.04;  0.89]       1.60 [1.04; 2.44]  ← beneficial
  w        -0.54 [-0.90; -0.17]       0.58 [0.41; 0.84]  ← harmful/control
  tau       0.59 [ 0.42;  0.82]
```

A 95% CrI that does **not** include 1 (for OR) or 0 (for log-OR) indicates
a credible effect.

### Interaction results (from `plot_interactions()` / `summarise_interactions()`)

| Statistic | Interpretation |
|-----------|----------------|
| `inclusion_prob` | **Posterior probability that the interaction is non-negligible.** Values ≥ 0.5 are conventionally "selected". Values near 0 mean the data strongly favour pure additivity for that pair. |
| `median_gamma` | Direction and magnitude of the interaction when included (log-OR scale). Negative = the combination is **less effective** than the sum of its parts (antagonism). Positive = **more effective** (synergy). |
| `OR_interaction` | exp(γ̂): the multiplicative modification of efficacy when both components co-occur, beyond their additive contributions. |

**Decision rule:** An interaction is considered meaningful if
`inclusion_prob ≥ 0.5` *and* the `OR_interaction` is clinically relevant
(e.g., > 1.2 or < 0.8 as a rough guide).

### Comparing models (additive vs. LASSO)

- If the LASSO and additive component ORs are **similar**, interactions add
  little and the simpler additive model may suffice.
- If they **diverge**, one or more interactions are being absorbed into the
  main effects in the additive model, and the LASSO model should be
  preferred for combination comparisons.

### Heterogeneity (τ)

`tau` is the between-study standard deviation of true log-ORs. Values
< 0.3 are generally considered small, 0.3–0.8 moderate, > 0.8 large (in
the context of psychological interventions; Turner et al. 2015).

---

## Accessing raw MCMC samples

```r
# Combined posterior samples (all chains stacked)
fit$samps_additive   # columns: d[1],...,d[Nc], ORd[1],...,ORd[Nc], tau
fit$samps_lasso      # columns: d[k], ORd[k], gamma[k], Ind[k], tau, eta

# Interaction summaries
fit$interaction_inclusion_prob   # named vector: P(Ind[k] = 1)
fit$interaction_effects          # named vector: median(gamma[k]), log-OR scale

# Custom posterior calculation
d_cr <- fit$samps_lasso[["d[3]"]]   # cr is component 3
cat("P(cr is beneficial):", mean(exp(d_cr) > 1), "\n")

# Posterior of a specific interaction
k <- place("cr", "th", Nc = length(components), components = components)
gamma_crth <- fit$samps_lasso[[paste0("gamma[", k, "]")]]
quantile(exp(gamma_crth), c(0.025, 0.5, 0.975))
```

---

## Interaction index utilities

```r
# Index of the (sr, sc) interaction
place("sr", "sc", Nc = length(components), components = components)

# Component names for interaction index 32
which.place(32, components)

# All interaction indices for a combination
all.interactions(c("se", "sd", "cr", "th"), Nc = length(components),
                 components = components)
```

---

## How to cite

If you use `bcnma` in a publication, please cite both the package and the
underlying methods paper:

```
Furukawa Y (2026). bcnma: Bayesian Component Network Meta-Analysis with
  LASSO Interactions. R package version 0.1.0.
  https://github.com/ykfrkw/bcnma

Furukawa Y, Sakata M, Yamamoto R, et al. (2024). Components and delivery
  formats of cognitive behavioral therapy for chronic insomnia in adults.
  JAMA Psychiatry, 81(2), 130-140.
  https://doi.org/10.1001/jamapsychiatry.2023.5060
```

---

## Acknowledgements

The SSVS spike-and-slab approach for component NMA was conceptually
developed with guidance from **Orestis Efthimiou** (University of Bern).
The statistical model follows the framework described in the CBTICNMA
analysis scripts.

---

## References

- Turner RM, et al. (2015). Predictive distributions for between-study
  heterogeneity and simple methods for their application in Bayesian
  meta-analysis. *Statistics in Medicine*, 34, 984–998.
- Higgins JPT, Whitehead A. (1996). Borrowing strength from external trials
  in a meta-analysis. *Statistics in Medicine*, 15, 2733–2749.
- George EI, McCulloch RE. (1993). Variable selection via Gibbs sampling.
  *Journal of the American Statistical Association*, 88, 881–889.
- Rücker G, Schwarzer G, et al. (2020). netmeta: Network meta-analysis
  using frequentist methods. R package.

---

## File structure

```
bcnma/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
│   ├── bcnma-package.R   ← package documentation
│   ├── generics.R        ← forest() generic
│   ├── bcnma.R           ← main bcnma() function
│   ├── helpers.R         ← data prep; place(), which.place(), all.interactions()
│   ├── jags_models.R     ← dynamic JAGS model string builders (internal)
│   └── methods.R         ← print, summary, forest, compare_combinations,
│                            plot_interactions, summarise_interactions
└── inst/extdata/
    ├── data_CBTICNMA.csv       ← bundled example dataset
    └── example_CBTICNMA.R     ← full worked example
```
