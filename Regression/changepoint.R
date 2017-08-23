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

m <- readRDS("Regression/changepoint.rds")
data <- list(N=100, D=c(rnorm(50, 100, 50), rnorm(50, 200, 50)))
f <- sampling(m, data=data, iter = 1000)
farr <- as.array(f)
mcmc_areas(farr, regex_pars = "^mu\\d")
mcmc_areas(farr, regex_pars = "^sigma\\d")
mcmc_dens_overlay(farr, regex_pars = "tau")

save(m2, file = "Regression/changepoint.RData")
