#' Perform a probabilistic PCA with ARD
#'
#' This function performs a probabilistic principal component analysis with
#' automatic relevance determination.
#'
#' @param xdf the data set to perform the analysis on
#' @param k the number of latent factors to consider
#' @param chains the number of chains to use which defaults to 2
#' @param iter the number of samples to pull which defaults to 1000
#' @param cores the number of cores to use which defaults to max(parallel::detectCores()-1, 1)
#' @param ... other arguments passed to the sampling method in rstan
#' @importFrom rstan sampling
#' @importFrom parallel detectCores
#' @return the stanfit from the sampled changepoint model
#' @export
#'
#' @examples
ppcaard <- function(xdf, k, chains=2, iter=500, cores=max(parallel::detectCores()-1, 1), ...)
{
  if(!inherits(xdf, "data.frame")) stop("Error: You need to supply a data.frame as xdf!")
  if(any(missing(k), length(k)>1)) stop("Error: You need to supply a scalar value for k!")
  data <- list(N=nrow(xdf), D=ncol(xdf), M=k)
  sampling(stanmodels$ppcaard, data=data, iter = iter, chains=chains, cores=cores, ...)
}
