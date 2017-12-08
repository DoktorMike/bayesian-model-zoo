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
#' \dontrun{
#'   data<-list(N=nrow(mtcars), y=mtcars$gear, X=mtcars[, c('drat', 'am')])
#'   sfit<-negbinom(data$y, data$X, iter = 500, chains=2, cores=2, allpars = TRUE)
#'   library(bayesplot)
#'   mcmc_combo(as.array(sfit), regex_pars = "tau")
#'   yrep<-as.matrix(sfit)
#'   ppc_dens_overlay(data$y, yrep[1:50,grep("y_", colnames(yrep))])
#' }
negbinom <- function(y, x, chains=2, iter=1000, allpars=FALSE, cores=max(parallel::detectCores()-1, 1), ...)
{
  data <- list(N=nrow(x), K=ncol(x), y=y, X=x)
  if(allpars) pars <- NA else pars <- c("beta", "mu", "phi", "y_rep")
  rstan::sampling(stanmodels$multiple_negative_binomial_regression, data=data, iter = iter, chains=chains, pars=pars, cores=cores, ...)
}
