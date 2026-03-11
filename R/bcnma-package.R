#' bcnma: Bayesian Component Network Meta-Analysis with LASSO Interactions
#'
#' Fits Bayesian random-effects component NMA for binary outcomes using JAGS.
#' Three interaction models are available: additive, LASSO (equiprobable SSVS),
#' and LASSO with informative priors. The interface mirrors \code{netmeta}/
#' \code{discomb} from the \pkg{meta} ecosystem.
#'
#' @section Main function:
#' \code{\link{bcnma}} — fit the model
#'
#' @section Output functions:
#' \code{\link{print.bcnma}}, \code{\link{summary.bcnma}},
#' \code{\link{forest.bcnma}}, \code{\link{plot_interactions}},
#' \code{\link{summarise_interactions}}, \code{\link{compare_combinations}}
#'
#' @section Index helpers:
#' \code{\link{place}}, \code{\link{which.place}}, \code{\link{all.interactions}}
#'
#' @references
#' Furukawa Y, et al. (2024). Components and delivery formats of cognitive
#' behavioral therapy for chronic insomnia in adults: A systematic review and
#' component network meta-analysis. \emph{JAMA Psychiatry}, 81(2), 130-140.
#' \doi{10.1001/jamapsychiatry.2023.5060}
#'
#' Turner RM, et al. (2015). Predictive distributions for between-study
#' heterogeneity. \emph{Statistics in Medicine}, 34, 984-998.
#'
#' @keywords internal
"_PACKAGE"
