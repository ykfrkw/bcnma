# =============================================================================
# bcnma — Worked example using CBTICNMA data
# Data source: Furukawa et al. (2024), JAMA Psychiatry
#              https://github.com/ykfrkw/CBTICNMA
# =============================================================================

library(bcnma)
library(ggplot2)

# ---- 1. Load data -----------------------------------------------------------
dat_raw <- read.csv(system.file("extdata", "data_CBTICNMA.csv", package = "bcnma_bin"))
names(dat_raw)[1] <- "study"  # fix BOM character in first column name

# Prepare arm-level data: one row per arm, aggregate duplicates if any
library(dplyr)
dat <- dat_raw |>
  dplyr::select(study, treatment_component, r, n) |>
  dplyr::group_by(study, treatment_component) |>
  dplyr::summarise(r = sum(r), n = sum(n), .groups = "drop")

# 17 components used in the original analysis
components <- c("se", "sd", "cr", "th", "cw", "sr", "sc", "re", "pi",
                "w",  "ns", "he", "tg", "ind", "gp", "ff", "ae")

# ---- 2a. Additive model (no interactions) -----------------------------------
fit_add <- bcnma(
  data         = dat,
  studlab      = study,
  treat        = treatment_component,
  event        = r,
  n            = n,
  components   = components,
  interactions = "none",
  n.chains     = 3,
  n.adapt      = 2000,
  n.burnin     = 5000,
  n.iter       = 20000,
  seed         = 42
)
print(fit_add)

# ---- 2b. LASSO model (all pairwise interactions, equiprobable) --------------
fit_lasso <- bcnma(
  data         = dat,
  studlab      = study,
  treat        = treatment_component,
  event        = r,
  n            = n,
  components   = components,
  interactions = "lasso",
  lasso_prob   = 0.5,
  g            = 100,
  n.chains     = 4,
  n.adapt      = 2000,
  n.burnin     = 5000,
  n.iter       = 20000,
  seed         = 42
)
print(fit_lasso)

# ---- 2c. LASSO informative: only clinically plausible interactions ----------
active <- c("se", "sd", "cr", "th", "sr", "sc", "re", "ff")

fit_inform <- bcnma(
  data              = dat,
  studlab           = study,
  treat             = treatment_component,
  event             = r,
  n                 = n,
  components        = components,
  interactions      = "lasso_informative",
  active_components = active,
  lasso_prob        = 0.8,
  g                 = 100,
  n.chains          = 4,
  n.adapt           = 2000,
  n.burnin          = 5000,
  n.iter            = 20000,
  seed              = 42
)
print(fit_inform)

# ---- 3. Visualise -----------------------------------------------------------

# Forest plot: additive vs. LASSO component OR estimates
forest(fit_lasso, model = "both")

# Interaction overview (gamma vs. inclusion probability)
plot_interactions(fit_lasso, top_n = 10)

# Full interaction table (sorted by |effect|)
inter_tbl <- summarise_interactions(fit_lasso, threshold = 0.5)
print(head(inter_tbl, 15))

# ---- 4. Compare combinations -----------------------------------------------

# Full CBT-I package vs. sleep hygiene only (se)
compare_combinations(fit_lasso,
  combo1 = c("se", "sd", "cr", "th", "sr", "sc"),
  combo2 = c("se"))

# Full package vs. no-component baseline
compare_combinations(fit_lasso,
  combo1 = c("se", "sd", "cr", "th", "sr", "sc"),
  combo2 = character(0))

# Same comparison, additive model only (no interaction adjustment)
compare_combinations(fit_lasso,
  combo1 = c("se", "sd", "cr", "th", "sr", "sc"),
  combo2 = c("se"),
  model  = "additive")

# ---- 5. Raw MCMC access -----------------------------------------------------

# Posterior OR for sleep restriction (sr = component 6)
d_sr <- fit_lasso$samps_lasso[["d[6]"]]
cat("SR vs. no component — OR [95% CrI]:",
    round(quantile(exp(d_sr), c(0.5, 0.025, 0.975)), 2), "\n")

# P(SR is beneficial, i.e. OR > 1)
cat("P(SR is beneficial):", round(mean(exp(d_sr) > 1), 3), "\n")

# Which interaction index corresponds to (sd, th)?
k_sdth <- place("sd", "th", Nc = length(components), components = components)
gamma_sdth <- fit_lasso$samps_lasso[[paste0("gamma[", k_sdth, "]")]]
cat("sd:th interaction OR [95% CrI]:",
    round(quantile(exp(gamma_sdth), c(0.5, 0.025, 0.975)), 3), "\n")
