# © 2026 Yuki Furukawa. All rights reserved.
# Author: Yuki Furukawa
# Date: 2026-03-11
#
# File description:
# S3 methods for bcnma_cont objects (continuous outcomes, SMD/MD scale).
# Mirrors the structure of methods.R for bcnma (binary outcomes).

# ---- print ----

#' @export
print.bcnma_cont <- function(x, digits = 2,
                              model = c("both", "additive", "lasso"), ...) {
  model <- match.arg(model)

  cat("Bayesian Component NMA — Continuous outcomes (bcnma_cont)\n")
  cat(rep("-", 60), "\n", sep = "")
  cat("Components    :", x$n_components, "\n")
  cat("Studies       :", x$Ns, "\n")
  cat("Effect measure:", x$sm, "\n")
  cat("Interactions  :", x$interactions_model, "\n")
  if (!is.null(x$reference.group)) cat("Reference     :", x$reference.group, "\n")
  cat("\n")

  show_add  <- model %in% c("both", "additive") && !is.null(x$summary_additive)
  show_lass <- model %in% c("both", "lasso")    && !is.null(x$summary_lasso)

  if (show_add) {
    cat("Component effects — Additive model\n")
    cat(rep("-", 60), "\n", sep = "")
    .print_cont_table(x$summary_additive, x$components, x$sm, digits)
  }

  if (show_lass) {
    cat(sprintf("\nComponent effects — LASSO model (%s)\n", x$interactions_model))
    cat(rep("-", 60), "\n", sep = "")
    .print_cont_table(x$summary_lasso, x$components, x$sm, digits)

    n_selected <- sum(x$interaction_inclusion_prob >= 0.5, na.rm = TRUE)
    cat(sprintf("\nInteractions with inclusion prob >= 0.5: %d / %d\n",
                n_selected, x$Ninter))
    if (n_selected > 0) {
      top <- summarise_interactions(x, threshold = 0.5, top_n = 10)
      top_sel <- top[top$selected, ]
      if (nrow(top_sel) > 0) {
        cat(sprintf("  %-20s  %8s  %8s\n", "Pair",
                    paste0("median ", x$sm), "P(incl)"))
        for (i in seq_len(nrow(top_sel))) {
          cat(sprintf("  %-20s  %8.3f  %8.3f\n",
                      top_sel$component_pair[i],
                      top_sel$median_gamma[i],
                      top_sel$inclusion_prob[i]))
        }
      }
    }
  }

  invisible(x)
}

# Internal: print component table for continuous models (no ORd column)
.print_cont_table <- function(smry, components, sm, digits) {
  Nc     <- length(components)
  d_rows <- smry[paste0("d[", seq_len(Nc), "]"), , drop = FALSE]
  tau_row <- smry["tau", , drop = FALSE]

  fmt <- function(m, lo, hi) sprintf("%.*f [%.*f; %.*f]",
                                      digits, m, digits, lo, digits, hi)

  cat(sprintf("  %-8s  %-28s\n", "Comp", paste0(sm, " [95% CrI]")))
  cat("  ", rep("-", 38), "\n", sep = "")
  for (k in seq_len(Nc)) {
    cat(sprintf("  %-8s  %-28s\n",
                components[k],
                fmt(d_rows[k, "50%"], d_rows[k, "2.5%"], d_rows[k, "97.5%"])))
  }
  if (!is.null(tau_row) && nrow(tau_row) > 0) {
    cat(sprintf("  %-8s  %-28s\n", "tau",
                fmt(tau_row[1, "50%"], tau_row[1, "2.5%"], tau_row[1, "97.5%"])))
  }
}

# ---- summary ----

#' @export
summary.bcnma_cont <- function(object, digits = 3, ...) {
  cat("=== Bayesian cNMA Summary (Continuous Outcomes) ===\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  cat("Data:\n")
  cat("  Studies      :", object$Ns, "\n")
  cat("  Components   :", object$n_components, "\n")
  cat("  Effect measure:", object$sm, "\n")
  cat("  Interactions :", object$Ninter, "pairs total\n")
  cat("  MCMC chains  :", object$mcmc_settings$n.chains,
      "  iter:", object$mcmc_settings$n.iter, "\n\n")

  print(object, digits = digits, model = "both")
  invisible(object)
}

# ---- forest ----

#' Forest plot of component SMD (or MD) effects
#'
#' @param x      bcnma_cont object
#' @param model  character: \code{"additive"}, \code{"lasso"}, or \code{"both"}
#' @param ref_line numeric: reference line (default 0 for SMD/MD)
#' @param xlab   character: x-axis label
#' @param title  character: plot title
#' @param ...    additional arguments (unused)
#' @return ggplot object (invisible)
#' @export
forest.bcnma_cont <- function(x,
                               model    = c("both", "additive", "lasso"),
                               ref_line = 0,
                               xlab     = NULL,
                               title    = "Component effects",
                               ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  model <- match.arg(model)
  if (is.null(xlab)) xlab <- paste0(x$sm, " vs. no component")

  Nc    <- x$n_components
  comps <- x$components

  extract_df <- function(smry, label) {
    d_rows <- smry[paste0("d[", seq_len(Nc), "]"), , drop = FALSE]
    data.frame(
      component = factor(comps, levels = rev(comps)),
      median    = d_rows[, "50%"],
      lo        = d_rows[, "2.5%"],
      hi        = d_rows[, "97.5%"],
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
    ggplot2::labs(x = xlab, y = NULL, title = title,
                  colour = "Model", shape = "Model") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}

# ---- compare_combinations for bcnma_cont ----

#' Compare posterior treatment effects for two component combinations
#' (continuous outcomes)
#'
#' Computes the posterior distribution of the SMD (or MD) for combination A
#' vs. combination B, accounting for main effects and (optionally) interaction
#' terms from the LASSO model.
#'
#' @param x            bcnma_cont object
#' @param combo1       character or integer vector: components in combination A
#' @param combo2       character or integer vector: components in combination B.
#'                     Use \code{character(0)} to compare against baseline.
#' @param model        character: \code{"lasso"} (default) or \code{"additive"}
#' @param probs        numeric vector: quantiles to return
#'
#' @return named list (same structure as \code{\link{compare_combinations}})
#' @export
compare_combinations.bcnma_cont <- function(x,
                                             combo1,
                                             combo2 = character(0),
                                             model  = c("lasso", "additive"),
                                             probs  = c(0.025, 0.5, 0.975)) {
  model <- match.arg(model)
  if (!inherits(x, "bcnma_cont")) stop("'x' must be a bcnma_cont object.")

  comps <- x$components
  Nc    <- x$n_components

  .resolve <- function(combo) {
    if (length(combo) == 0) return(integer(0))
    if (is.character(combo)) {
      unknown <- setdiff(combo, comps)
      if (length(unknown) > 0) {
        stop("Unknown components: ", paste(unknown, collapse = ", "))
      }
      which(comps %in% combo)
    } else {
      idx <- as.integer(combo)
      if (any(idx < 1 | idx > Nc)) stop("Indices out of range [1, ", Nc, "]")
      idx
    }
  }

  idx1 <- .resolve(combo1)
  idx2 <- .resolve(combo2)
  if (length(idx1) == 0) stop("'combo1' must contain at least one component.")

  if (model == "lasso") {
    if (is.null(x$samps_lasso)) {
      warning("LASSO samples not available; falling back to additive model.")
      model <- "additive"
      samps <- x$samps_additive
    } else {
      samps <- x$samps_lasso
    }
  } else {
    samps <- x$samps_additive
  }
  if (is.null(samps)) stop("No MCMC samples available for model = '", model, "'.")

  n_samp <- nrow(samps)
  effect <- numeric(n_samp)

  for (j in idx1) {
    col <- paste0("d[", j, "]")
    if (!col %in% colnames(samps)) stop("Column '", col, "' not found.")
    effect <- effect + samps[[col]]
  }
  for (j in idx2) {
    col <- paste0("d[", j, "]")
    if (!col %in% colnames(samps)) stop("Column '", col, "' not found.")
    effect <- effect - samps[[col]]
  }

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

  combo1_names <- if (is.character(combo1) && length(combo1) > 0) combo1 else comps[idx1]
  combo2_names <- if (length(combo2) == 0) "baseline" else
    (if (is.character(combo2) && length(combo2) > 0) combo2 else comps[idx2])

  cat(sprintf("Comparison: [%s] vs. [%s] | model: %s\n",
              paste(combo1_names, collapse = "+"),
              if (is.character(combo2_names)) paste(combo2_names, collapse = "+") else combo2_names,
              model))
  est <- quantile(effect, probs)
  cat(sprintf("  %s: %.3f [%.3f; %.3f]\n", x$sm, est[2], est[1], est[3]))

  invisible(list(
    estimate       = est,
    combo1         = combo1_names,
    combo2         = combo2_names,
    model          = model,
    sm             = x$sm,
    n_mcmc         = n_samp,
    effect_samples = effect
  ))
}
