# © 2026 Yuki Furukawa. All rights reserved.
# Author: Yuki Furukawa
# Date: 2026-03-09
# Contact: furukawa.y.psy@gmail.com
#
# File description:
# Helper functions for Bayesian component NMA with LASSO/SSVS interactions.
# Includes arm-level data preparation and interaction indexing utilities.

# ---- Data preparation ----

#' Prepare arm-level binary data for Bayesian cNMA
#'
#' Converts long-format arm-level data into JAGS-ready matrices.
#'
#' @param data        data.frame with study-arm level rows
#' @param studlab_col character: column name for study ID
#' @param treat_col   character: column name for treatment combination string
#' @param event_col   character: column name for event count
#' @param n_col       character: column name for sample size
#' @param components  character vector of component abbreviations (order matters)
#' @param component_separator character: separator used in treat column (default "+")
#' @return list with JAGS data and metadata
.bcnma_prepare_data <- function(data, studlab_col, treat_col, event_col, n_col,
                                components, component_separator = "+") {
  d <- data.frame(
    study = as.character(data[[studlab_col]]),
    treat = as.character(data[[treat_col]]),
    r     = as.integer(floor(as.numeric(data[[event_col]]))),
    n     = as.integer(floor(as.numeric(data[[n_col]])))
  )

  d <- d[order(d$study), ]
  studies  <- unique(d$study)
  Ns       <- length(studies)
  Nc       <- length(components)

  # Assign arm numbers within each study
  na_vec <- integer(Ns)
  d$arm  <- 0L
  for (i in seq_along(studies)) {
    idx       <- which(d$study == studies[i])
    na_vec[i] <- length(idx)
    d$arm[idx] <- seq_along(idx)
  }

  max_arms <- max(na_vec)

  # Component indicator array: xA[study, arm, component]
  # xA[i, k, j] = 1 if component j is present in arm k of study i
  xA <- array(0L, dim = c(Ns, max_arms, Nc))
  for (i in seq_along(studies)) {
    rows <- d[d$study == studies[i], ]
    for (k in seq_len(na_vec[i])) {
      parts <- trimws(strsplit(rows$treat[k], component_separator, fixed = TRUE)[[1]])
      for (j in seq_along(components)) {
        xA[i, k, j] <- as.integer(components[j] %in% parts)
      }
    }
  }

  # Reference arm (arm 1) indicators: xB[study, component]
  xB <- xA[, 1, , drop = FALSE]
  dim(xB) <- c(Ns, Nc)

  # Reshape r and n into (Ns x max_arms), NA for unused arms
  r_mat <- matrix(NA_integer_, Ns, max_arms)
  n_mat <- matrix(NA_integer_, Ns, max_arms)
  for (i in seq_along(studies)) {
    rows <- d[d$study == studies[i], ]
    r_mat[i, seq_len(na_vec[i])] <- rows$r
    n_mat[i, seq_len(na_vec[i])] <- rows$n
  }

  # Check for components not found in any arm
  comp_used <- apply(xA, 3, sum)
  if (any(comp_used == 0)) {
    warning("The following components appear in no arm: ",
            paste(components[comp_used == 0], collapse = ", "))
  }

  list(
    r          = r_mat,
    n          = n_mat,
    xA         = xA,
    xB         = xB,
    na         = na_vec,
    Ns         = Ns,
    Nc         = Nc,
    max_arms   = max_arms,
    studies    = studies,
    components = components,
    data       = d
  )
}

#' Build pairwise interaction indicator array for JAGS
#'
#' @param xA array [Ns, max_arms, Nc] binary component indicators
#' @param Ns integer number of studies
#' @param max_arms integer maximum arms per study
#' @param Nc integer number of components
#' @return array [Ns, max_arms, Ninter] where Ninter = Nc*(Nc-1)/2
.build_interaction_array <- function(xA, Ns, max_arms, Nc) {
  Ninter <- Nc * (Nc - 1L) / 2L
  inter  <- array(0L, dim = c(Ns, max_arms, Ninter))
  idx    <- 0L
  for (p in seq_len(Nc - 1L)) {
    for (q in (p + 1L):Nc) {
      idx            <- idx + 1L
      inter[, , idx] <- xA[, , p] * xA[, , q]
    }
  }
  inter
}

# ---- Interaction indexing ----

#' Map component pair (i, j) to pairwise interaction index
#'
#' Numbering follows upper-triangle row-major order.
#' e.g., for Nc=4: (1,2)->1, (1,3)->2, (1,4)->3, (2,3)->4, (2,4)->5, (3,4)->6
#'
#' @param i          integer or character: first component (index or name)
#' @param j          integer or character: second component (index or name)
#' @param Nc         integer total number of components
#' @param components character vector of component names (required when i or j is character)
#' @return integer interaction index
place <- function(i, j, Nc, components = NULL) {
  .resolve_idx <- function(x) {
    if (is.character(x)) {
      if (is.null(components)) stop("'components' must be provided when i or j is character")
      idx <- which(components == x)
      if (length(idx) == 0) stop("Component '", x, "' not found in components vector")
      idx
    } else {
      as.integer(x)
    }
  }
  i <- .resolve_idx(i)
  j <- .resolve_idx(j)
  if (i == j) stop("i and j must refer to different components")
  if (i > j) { tmp <- i; i <- j; j <- tmp }
  as.integer((2L * Nc - i) * (i - 1L) / 2L + j - i)
}

#' Identify component pair from an interaction index
#'
#' @param idx        integer interaction index
#' @param components character vector of component names
#' @return named character vector: c(indices = "i j", names = "comp_i - comp_j")
which.place <- function(idx, components) {
  Nc <- length(components)
  for (k in seq_len(Nc - 1L)) {
    for (l in (k + 1L):Nc) {
      if (place(k, l, Nc) == idx) {
        return(c(
          indices = paste(k, l),
          names   = paste(components[k], "-", components[l])
        ))
      }
    }
  }
  stop(sprintf("Interaction index %d not found for Nc = %d", idx, Nc))
}

#' Get all pairwise interaction indices for a component combination
#'
#' @param combination integer or character vector of component indices (or names)
#' @param Nc          integer total number of components
#' @param components  character vector of component names (required if combination is character)
#' @return integer vector of interaction indices; empty if fewer than 2 components
all.interactions <- function(combination, Nc, components = NULL) {
  # Accept component names
  if (is.character(combination)) {
    if (is.null(components)) stop("'components' must be provided when 'combination' is character")
    combination <- which(components %in% combination)
  }
  combination <- sort(unique(as.integer(combination)))
  nc <- length(combination)
  if (nc < 2L) return(integer(0))
  out <- vector("integer", nc * (nc - 1L) / 2L)
  idx <- 0L
  for (i in seq_len(nc - 1L)) {
    for (j in (i + 1L):nc) {
      idx      <- idx + 1L
      out[idx] <- place(combination[i], combination[j], Nc)
    }
  }
  out
}

#' Summarise interaction results from a bcnma object
#'
#' @param x          bcnma object
#' @param threshold  numeric: inclusion probability threshold for "selected" label (default 0.5)
#' @param top_n      integer: how many top interactions to return (default all)
#' @return data.frame sorted by |median_gamma| descending
summarise_interactions <- function(x, threshold = 0.5, top_n = NULL) {
  if (is.null(x$interaction_inclusion_prob)) {
    stop("No interaction results found. Run bcnma() with interactions = 'lasso' or 'lasso_informative'.")
  }

  Nc     <- x$n_components
  Ninter <- Nc * (Nc - 1L) / 2L
  comps  <- x$components

  rows <- lapply(seq_len(Ninter), function(k) {
    pair <- which.place(k, comps)
    data.frame(
      index          = k,
      component_pair = pair[["names"]],
      median_gamma   = x$interaction_effects[k],
      OR_interaction = exp(x$interaction_effects[k]),
      inclusion_prob = x$interaction_inclusion_prob[k],
      selected       = x$interaction_inclusion_prob[k] >= threshold,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(abs(out$median_gamma), decreasing = TRUE), ]
  rownames(out) <- NULL

  if (!is.null(top_n)) out <- head(out, top_n)
  out
}
