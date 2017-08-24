library(HDInterval)
library(rstan)
library(bayesplot)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

hasChangepoint <- function(x, chains=2, iter=1000)
{
  m <- readRDS("Regression/changepoint.rds")
  data <- list(N=length(x), y=x)
  f <- sampling(m, data=data, iter = iter, chains=chains)
  fdf <- as.data.frame(f)
  resdf <- tibble(param=c("mu", "sigma"),
                  overlap=c(distOverlap(fdf$mu1, fdf$mu2), distOverlap(fdf$sigma1, fdf$sigma2)),
                  probability=1-overlap,
                  changepoint=median(fdf$tau))
  return(list(params=fdf, results=resdf))
}

m <- readRDS("Regression/changepoint.rds")
data <- list(N=100, y=c(rnorm(50, 100, 100), rnorm(50, 150, 100)))
f <- sampling(m, data=data, iter = 1000, chains=2)
# f <- vb(m, data=data, iter = 1000, algorithm="meanfield")
farr <- as.array(f)
mcmc_areas(farr, regex_pars = "^mu\\d")
mcmc_areas(farr, regex_pars = "^sigma\\d")
mcmc_dens_overlay(farr, regex_pars = "tau")
