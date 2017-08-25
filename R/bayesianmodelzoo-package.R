#' Applied Regression Modeling via RStan
#'
#' @docType package
#' @name bayesianmodelzoo-package
#' @aliases bayesmodelzoo
#' @useDynLib rstanarm, .registration = TRUE
#'
#' @import methods
#' @importFrom rstan optimizing sampling vb constrain_pars extract
#'   extract_sparse_parts get_posterior_mean stanc
#' @import stats
#' @import Rcpp
#' @import bayesplot
#' @import shinystan
#' @import rstantools
#' @export log_lik posterior_linpred posterior_predict posterior_interval
#' @export predictive_interval predictive_error prior_summary bayes_R2
#' @export loo_linpred loo_predict loo_predictive_interval
#' @export loo waic
#' @export launch_shinystan
#'
#' @description
#' \if{html}{
#'    \figure{stanlogo.png}{options: width="50px" alt="mc-stan.org"}
#'    \emph{Stan Development Team}
#' }
#'
#' The \pkg{bayesmodelzoo} package is an appendage to the \pkg{rstan} package
#'
NULL
