#' Quantify distribution overlap
#'
#' Measures the amount of overlap between two sampled numerical distributions x1
#' and x2. The return value states the degree of overlap. If 1 then 100 percent
#' is overlapped meanwhile a 0 indicates no overlap.
#'
#' @param x1 samples from distribution 1
#' @param x2 samples from distribution 2
#' @param breaks manually specify the breaks to use which if set to NULL (the default) will be estimated by histogram
#' @importFrom HDInterval hdi
#' @importFrom graphics hist
#' @return the overlap in percent represented as a number between 0 and 1
#' @export
#'
#' @examples
#' distOverlap(rnorm(10, 10, 5), rnorm(10, 5, 5))
distOverlap <- function(x1, x2, breaks=NULL) {
  if(!all(is.numeric(x2), is.numeric(x1))) stop("Both vectors supplied must be numeric!")
  if(is.null(breaks)) {
    nbreaks <- ceiling(diff(range(x1)) / as.numeric(diff(hdi(x1))/18))
    breaks <- seq(from=min(x1), to=max(x1), length.out=nbreaks)
  }
  # histograms
  x1Info <- hist(x1, plot = FALSE, breaks=breaks)
  x2Info <- hist(x2, breaks=c(-Inf, breaks, Inf), plot = FALSE)$density[2:length(breaks)]
  # get the overlap
  minHt <- pmin(x2Info, x1Info$density)
  overlap <- sum(minHt * diff(x1Info$breaks))

  return(overlap)
}

#' Investigate the possibility of a changepoint
#'
#' Checks whether it's probable that there is a changepoint in the unidimensional
#' time series x. Currently only one changepoint is supported. The model works by
#' investigating the hypothesis that there is a point in the time series at which
#' the generating probability distribution shifts. These distributions are assumed
#' to be gaussian with a mu1 and mu2. If mu1 and mu2 overlap completely there is no
#' changepoint. The distribution of the location of the changepoint is in the "tau"
#' parameter and is returned in the results as median.
#'
#' @param x the numeric vector representing the timeseries to investigate
#' @param chains number of MCMC chains to run defualts to 2
#' @param iter number of samples to pull from MCMC defaults to 1000
#' @param allpars return all parameter samples from the changepoint model which
#' defaults to FALSE and thus returns only the relevant ones for the decision
#' @param ... other arguments you want to send to sampling
#' @importFrom tibble tibble
#' @importFrom rstan sampling
#' @importFrom stats median
#' @return a list of names results containing a tibble with overlap probabilities
#' for mu and sigma as well as the samples for the relevant parameters of the model.
#' @export
#'
#' @examples
#' a <- hasChangepoint(c(rnorm(50, 100, 100), rnorm(50, 150, 100)))
#' print(a$results)
hasChangepoint <- function(x, chains=2, iter=1000, allpars=FALSE)
{
  # m <- readRDS("Regression/changepoint.rds")
  data <- list(N=length(x), y=x)
  if(allpars) pars <- NA else pars <- c("mu1", "mu2", "sigma1", "sigma2", "tau")
  f <- sampling(m, data=data, iter = iter, chains=chains, pars=pars)
  fdf <- as.data.frame(f)
  resdf <- tibble(param=c("mu", "sigma"),
                  overlap=c(distOverlap(fdf$mu1, fdf$mu2), distOverlap(fdf$sigma1, fdf$sigma2)),
                  probability=1-overlap,
                  changepoint=median(fdf$tau))
  return(list(params=fdf, results=resdf))
}
