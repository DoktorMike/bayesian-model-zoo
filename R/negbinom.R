#' Bayesian negative binomial regression
#'
#' This implements a Bayesian negative binomial regression by using the 2nd
#' parameterization from the Stan manual.
#'
#' @param y the response variable to uses
#' @param x the matrix of covariates to use
#' @param chains the number of chains to use which defaults to 2
#' @param iter the number of samples to pull which defaults to 1000
#' @param allpars decides if all parameters should be returned or not. FALSE which is
#' the default only gives the mu's, beta's, phi and predicted y_rep
#' @param cores the number of cores to use which defaults to max(parallel::detectCores()-1, 1)
#' @param ... other arguments passed to the sampling method in rstan
#'
#' @return the stanfit from the sampled model
#' @export
#'
#' @examples
#' a<-1
negbinom <- function(y, x, chains=2, iter=1000, allpars=FALSE, cores=max(parallel::detectCores()-1, 1), ...)
{
  data <- list(N=length(x), y=y, x=x)
  if(allpars) pars <- NA else pars <- c("beta", "mu", "phi", "y_rep")
  sampling(stanmodels$multiple_negative_binomial_regression, data=data, iter = iter, chains=chains, pars=pars, cores=cores, ...)
}
