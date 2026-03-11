# © 2026 Yuki Furukawa. All rights reserved.
# Author: Yuki Furukawa
# Date: 2026-03-09
# Contact: furukawa.y.psy@gmail.com
#
# File description:
# S3 methods and standalone functions for bcnma objects:
#   print.bcnma       — component effects table
#   summary.bcnma     — full results including interactions
#   forest.bcnma      — forest plot of component OR
#   compare_combinations — posterior distribution for any combination pair

# ---- print ----

#' @export
print.bcnma <- function(x, digits = 2, model = c("both", "additive", "lasso"), ...) {
  model <- match.arg(model)

  cat("Bayesian Component NMA (bcnma)\n")
  cat(rep("-", 55), "\n", sep = "")
  cat("Components    :", x$n_components, "\n")
  cat("Studies       :", x$Ns, "\n")
  cat("Effect measure: OR (log-OR scale in d[k])\n")
  cat("Interactions  :", x$interactions_model, "\n")
  if (!is.null(x$reference.group)) cat("Reference     :", x$reference.group, "\n")
  cat("\n")

  show_add  <- model %in% c("both", "additive") && !is.null(x$summary_additive)
  show_lass <- model %in% c("both", "lasso")    && !is.null(x$summary_lasso)

  if (show_add) {
    cat("Component effects — Additive model\n")
    cat(rep("-", 55), "\n", sep = "")
    .print_component_table(x$summary_additive, x$components, digits)
  }

  if (show_lass) {
    cat("\nComponent effects — LASSO model (", x$interactions_model, ")\n", sep = "")
    cat(rep("-", 55), "\n", sep = "")
    .print_component_table(x$summary_lasso, x$components, digits)

    # Top selected interactions
    n_selected <- sum(x$interaction_inclusion_prob >= 0.5, na.rm = TRUE)
    cat(sprintf("\nInteractions with inclusion prob >= 0.5: %d / %d\n",
                n_selected, x$Ninter))
    if (n_selected > 0) {
      top <- summarise_interactions(x, threshold = 0.5, top_n = 10)
      top_sel <- top[top$selected, ]
      if (nrow(top_sel) > 0) {
        cat(sprintf("  %-20s  %6s  %8s  %8s\n",
                    "Pair", "median OR", "P(incl)", "selected"))
        for (i in seq_len(nrow(top_sel))) {
          cat(sprintf("  %-20s  %6.3f  %8.3f  %s\n",
                      top_sel$component_pair[i],
                      top_sel$OR_interaction[i],
                      top_sel$inclusion_prob[i],
                      if (top_sel$selected[i]) "*" else ""))
        }
      }
    }
  }

  invisible(x)
}

# Internal: print component effects table from MCMCsummary
.print_component_table <- function(smry, components, digits) {
  Nc <- length(components)
  # Find d[k] rows (log-OR)
  d_rows  <- smry[paste0("d[",   seq_len(Nc), "]"), , drop = FALSE]
  or_rows <- smry[paste0("ORd[", seq_len(Nc), "]"), , drop = FALSE]

  tau_row <- smry["tau", , drop = FALSE]

  fmt <- function(m, lo, hi) sprintf("%.*f [%.*f; %.*f]",
                                      digits, m, digits, lo, digits, hi)

  cat(sprintf("  %-8s  %-26s  %-26s\n", "Comp", "log-OR [95% CrI]", "OR [95% CrI]"))
  cat("  ", rep("-", 62), "\n", sep = "")
  for (k in seq_len(Nc)) {
    cat(sprintf("  %-8s  %-26s  %-26s\n",
                components[k],
                fmt(d_rows[k, "50%"], d_rows[k, "2.5%"], d_rows[k, "97.5%"]),
                fmt(or_rows[k, "50%"], or_rows[k, "2.5%"], or_rows[k, "97.5%"])
    ))
  }
  if (!is.null(tau_row) && nrow(tau_row) > 0) {
    cat(sprintf("  %-8s  %-26s\n", "tau",
                fmt(tau_row[1, "50%"], tau_row[1, "2.5%"], tau_row[1, "97.5%"])))
  }
}

# ---- summary ----

#' @export
summary.bcnma <- function(object, digits = 3, ...) {
  cat("=== Bayesian cNMA Summary ===\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Data:\n")
  cat("  Studies      :", object$Ns, "\n")
  cat("  Components   :", object$n_components, "\n")
  cat("  Interactions :", object$Ninter, "pairs total\n")
  cat("  MCMC chains  :", object$mcmc_settings$n.chains,
      "  iter:", object$mcmc_settings$n.iter,
      "  thin:", object$mcmc_settings$n.thin, "\n\n")

  print(object, digits = digits, model = "both")

  cat("\nInteraction model: ", object$interactions_model, "\n")
  if (object$interactions_model == "lasso_informative") {
    cat("Active components:", paste(object$active_components, collapse = ", "), "\n")
    cat("Modelled interactions:", length(object$included_interactions),
        "/ fixed to 0:", object$Ninter - length(object$included_interactions), "\n")
  }

  invisible(object)
}

# ---- forest ----

#' Forest plot of component OR effects
#'
#' @param x      bcnma object
#' @param model  character: `"additive"`, `"lasso"`, or `"both"` (default)
#' @param ref_line numeric: reference line position (default 1 for OR)
#' @param xlab   character: x-axis label
#' @param title  character: plot title
#' @param ...    additional arguments (unused)
#' @return ggplot object (invisible)
#' @export
forest.bcnma <- function(x,
                          model    = c("both", "additive", "lasso"),
                          ref_line = 1,
                          xlab     = "Odds Ratio vs. no component",
                          title    = "Component effects",
                          ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  model <- match.arg(model)

  Nc     <- x$n_components
  comps  <- x$components

  extract_df <- function(smry, label) {
    or_rows <- smry[paste0("ORd[", seq_len(Nc), "]"), , drop = FALSE]
    data.frame(
      component = factor(comps, levels = rev(comps)),
      median    = or_rows[, "50%"],
      lo        = or_rows[, "2.5%"],
      hi        = or_rows[, "97.5%"],
      model     = label,
      stringsAsFactors = FALSE
    )
  }

  dfs <- list()
  if (model %in% c("additive", "both") && !is.null(x$summary_additive)) {
    dfs[["Additive"]] <- extract_df(x$summary_additive, "Additive")
  }
  if (model %in% c("lasso", "both") && !is.null(x$summary_lasso)) {
    lbl <- if (x$interactions_model == "lasso") "LASSO" else "LASSO (informative)"
    dfs[[lbl]] <- extract_df(x$summary_lasso, lbl)
  }

  if (length(dfs) == 0) stop("No model results available to plot.")
  df <- do.call(rbind, dfs)

  dodge <- if (length(dfs) > 1) ggplot2::position_dodge(width = 0.5) else ggplot2::position_identity()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = median, y = component,
                                         xmin = lo, xmax = hi,
                                         colour = model, shape = model)) +
    ggplot2::geom_vline(xintercept = ref_line, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(height = 0.25), position = dodge) +
    ggplot2::geom_point(size = 3, position = dodge) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = xlab, y = NULL, title = title,
                  colour = "Model", shape = "Model") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}

# ---- compare_combinations ----

#' Compare posterior treatment effects for two component combinations
#'
#' Computes the posterior distribution of the OR (or log-OR) for combination A
#' vs combination B, accounting for main effects and (optionally) pairwise
#' interaction terms estimated by the LASSO/SSVS model.
#'
#' @param x            bcnma object
#' @param combo1       character vector of component names (or integer indices) for combination A
#' @param combo2       character vector of component names (or integer indices) for combination B.
#'                     Use `character(0)` or `integer(0)` to compare against "no components" (baseline).
#' @param model        character: `"lasso"` (default) or `"additive"`
#' @param exponentiate logical: return OR (`TRUE`, default) or log-OR (`FALSE`)
#' @param probs        numeric vector: quantiles to return (default 95% CrI + median)
#'
#' @return named list:
#' \describe{
#'   \item{estimate}{named numeric: posterior quantiles}
#'   \item{combo1}{character: combination A component names}
#'   \item{combo2}{character: combination B component names (or "baseline")}
#'   \item{model}{character: model used}
#'   \item{n_mcmc}{integer: number of posterior samples used}
#'   \item{effect_samples}{numeric vector: full posterior sample (for further use)}
#' }
#'
#' @examples
#' \dontrun{
#' # Full CBT-I package vs. sleep hygiene only
#' compare_combinations(fit,
#'   combo1 = c("se", "sd", "cr", "th", "sr", "sc"),
#'   combo2 = c("se"))
#'
#' # Single component vs. baseline
#' compare_combinations(fit, combo1 = "sr", combo2 = character(0))
#' }
#' @export
compare_combinations <- function(x,
                                  combo1,
                                  combo2       = character(0),
                                  model        = c("lasso", "additive"),
                                  exponentiate = TRUE,
                                  probs        = c(0.025, 0.5, 0.975)) {
  model <- match.arg(model)

  if (!inherits(x, "bcnma")) stop("'x' must be a bcnma object.")

  comps <- x$components
  Nc    <- x$n_components

  # Resolve component names → integer indices
  .resolve <- function(combo) {
    if (length(combo) == 0) return(integer(0))
    if (is.character(combo)) {
      unknown <- setdiff(combo, comps)
      if (length(unknown) > 0) {
        stop("Unknown components: ", paste(unknown, collapse = ", "),
             "\nAvailable: ", paste(comps, collapse = ", "))
      }
      which(comps %in% combo)
    } else {
      idx <- as.integer(combo)
      if (any(idx < 1 | idx > Nc)) stop("Component indices out of range [1, ", Nc, "]")
      idx
    }
  }

  idx1 <- .resolve(combo1)
  idx2 <- .resolve(combo2)

  if (length(idx1) == 0) stop("'combo1' must contain at least one component.")

  # Select MCMC sample data frame
  if (model == "lasso") {
    if (is.null(x$samps_lasso)) {
      warning("LASSO samples not available; falling back to additive model.")
      model  <- "additive"
      samps  <- x$samps_additive
    } else {
      samps  <- x$samps_lasso
    }
  } else {
    samps  <- x$samps_additive
  }
  if (is.null(samps)) stop("No MCMC samples available for model = '", model, "'.")

  n_samp <- nrow(samps)
  effect <- numeric(n_samp)

  # ---- Main effects ----
  for (j in idx1) {
    col <- paste0("d[", j, "]")
    if (!col %in% colnames(samps)) stop("Column '", col, "' not found in MCMC samples.")
    effect <- effect + samps[[col]]
  }
  for (j in idx2) {
    col <- paste0("d[", j, "]")
    if (!col %in% colnames(samps)) stop("Column '", col, "' not found in MCMC samples.")
    effect <- effect - samps[[col]]
  }

  # ---- Interaction effects (LASSO model only) ----
  if (model == "lasso" && !is.null(x$samps_lasso)) {
    inter1 <- all.interactions(idx1, Nc)
    inter2 <- all.interactions(idx2, Nc)

    for (k in inter1) {
      col <- paste0("gamma[", k, "]")
      if (col %in% colnames(samps)) effect <- effect + samps[[col]]
    }
    for (k in inter2) {
      col <- paste0("gamma[", k, "]")
      if (col %in% colnames(samps)) effect <- effect - samps[[col]]
    }
  }

  if (exponentiate) effect <- exp(effect)

  combo1_names <- if (is.character(combo1) && length(combo1) > 0) combo1 else comps[idx1]
  combo2_names <- if (length(combo2) == 0) "baseline" else
    (if (is.character(combo2) && length(combo2) > 0) combo2 else comps[idx2])

  cat(sprintf("Comparison: [%s] vs. [%s] | model: %s\n",
              paste(combo1_names, collapse = "+"),
              if (is.character(combo2_names)) paste(combo2_names, collapse = "+") else combo2_names,
              model))
  cat(sprintf("  %s: ", if (exponentiate) "OR" else "log-OR"))
  est <- quantile(effect, probs)
  cat(sprintf("%.3f [%.3f; %.3f]\n", est[2], est[1], est[3]))

  invisible(list(
    estimate       = est,
    combo1         = combo1_names,
    combo2         = combo2_names,
    model          = model,
    exponentiate   = exponentiate,
    n_mcmc         = n_samp,
    effect_samples = effect
  ))
}

#' Plot LASSO interaction overview: coefficient vs. inclusion probability
#'
#' Scatter plot of median gamma (x-axis) vs. mean P(Ind=1) (y-axis) for all
#' pairwise interactions, labelling the top interactions by |gamma|.
#'
#' @param x       bcnma object
#' @param top_n   integer: number of interaction pairs to label (default 10)
#' @param threshold numeric: draw horizontal reference line at this inclusion prob (default 0.5)
#' @return ggplot object (invisible)
#' @export
plot_interactions <- function(x, top_n = 10, threshold = 0.5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (is.null(x$interaction_effects)) {
    stop("No interaction results. Run bcnma() with interactions = 'lasso' or 'lasso_informative'.")
  }

  df <- summarise_interactions(x)

  top_labels <- head(df[order(abs(df$median_gamma), decreasing = TRUE), ], top_n)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = median_gamma, y = inclusion_prob,
                                         label = component_pair)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "red") +
    ggplot2::geom_point(ggplot2::aes(colour = selected), size = 2, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#d62728", "FALSE" = "grey60"),
                                  labels = c("TRUE" = paste0(">= ", threshold),
                                             "FALSE" = paste0("< ", threshold)),
                                  name = "Selected") +
    ggplot2::labs(
      x     = "Median interaction coefficient (log-OR scale)",
      y     = "Inclusion probability P(Ind = 1)",
      title = "LASSO interaction overview"
    ) +
    ggplot2::theme_bw(base_size = 12)

  # Add labels for top interactions (requires ggrepel if available)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = top_labels,
      ggplot2::aes(label = component_pair),
      size = 3, max.overlaps = 20
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = top_labels,
      ggplot2::aes(label = component_pair),
      size = 3, vjust = -0.8
    )
  }

  print(p)
  invisible(p)
}
