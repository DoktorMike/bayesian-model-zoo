
logsat <- function(x, a) log(x/a+1)
fisk <- function(x, a, b) 1/(1+(x/a)^(-b) )


load("data/bayesVsFreqAgain.RData")
write.csv(mydf, file = "data/data.csv")



# First model -------------------------------------------------------------
library(rstan)
library(bayesplot)
library(ggplot2)
library(dautility)
library(HDInterval)

m1 <- stan_model(file = "Regression/mmm_linear_1.stan")
data1 <- list(N=nrow(mydf),
              Kxmi=select(mydf, matches("impressions_")) %>% ncol(),
              Kxmc=select(mydf, clicks_Online_paidSearch_All) %>% ncol(),
              Kxmg=select(mydf, matches("gross_")) %>% ncol(),
              Kxwa=select(mydf, precipitation, wind_speed, mean_temperature) %>% ncol(),
              y=mydf$newusers,
              xswd=select(mydf, matches("WDay")),
              xsm=microSeason2(min(mydf$date), max(mydf$date), agglev = "daily", season = "yearly_monthly")[,-1],
              xmi=select(mydf, matches("impressions_")) %>% setNames(tolower(gsub("_Online", "", names(.)))),
              xmc=select(mydf, clicks_Online_paidSearch_All),
              xmg=select(mydf, matches("gross_")),
              xwa=select(mydf, precipitation, wind_speed, mean_temperature)
)
f1 <- sampling(m1, data=data1, iter = 500)
f1arr <- as.array(f1)
f1df <- as.data.frame(f1)
ppc_dens_overlay(y=data1$y, yrep = f1arr[1:23,1,grep("^yhat\\[*", names(f1), value = T)])
ppc_stat(y=data1$y, yrep = f1arr[1:23,1,grep("^yhat\\[*", names(f1), value = T)])
mcmc_areas(f1arr[,1,grep("^bmg\\[*", names(f1), value = T)])

tmpnames <- grep("^bmg\\[*", names(f1), value = T)
a<-lapply(1:length(tmpnames), function(x) data1$xmg[,x] %o% f1df[, tmpnames[x]])



# Changepoint test --------------------------------------------------------
library(HDInterval)
postPriorOverlap <-
  function( paramSampleVec, prior, ..., yaxt="n", ylab="",
            xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
            xlim=range(paramSampleVec), breaks=NULL) {

    # Does a posterior histogram for a single parameter, adds the prior,
    #   displays and calculates the overlap.
    # Returns the overlap.

    oldpar <- par(xpd=NA) ; on.exit(par(oldpar))

    # get breaks: a sensible number over the hdi; cover the full range (and no more);
    #   equal spacing.
    if (is.null(breaks)) {
      nbreaks <- ceiling(diff(range(paramSampleVec)) / as.numeric(diff(hdi(paramSampleVec))/18))
      breaks <- seq(from=min(paramSampleVec), to=max(paramSampleVec), length.out=nbreaks)
    }
    # plot posterior histogram.
    histinfo <- hist(paramSampleVec, xlab=xlab, yaxt=yaxt, ylab=ylab,
                     freq=FALSE, border='white', col='skyblue',
                     xlim=xlim, main=main, cex=cex, cex.lab=cex.lab,
                     breaks=breaks)

    if (is.numeric(prior))  {
      # plot the prior if it's numeric
      priorInfo <- hist(prior, breaks=c(-Inf, breaks, Inf), add=TRUE,
                        freq=FALSE, col='yellow', border='white')$density[2:length(breaks)]
    } else if (is.function(prior)) {
      if(class(try(prior(0.5, ...), TRUE)) == "try-error")
        stop(paste("Incorrect arguments for the density function", substitute(prior)))
      priorInfo <- prior(histinfo$mids, ...)
    }
    # get (and plot) the overlap
    minHt <- pmin(priorInfo, histinfo$density)
    rect(breaks[-length(breaks)], rep(0, length(breaks)-1), breaks[-1], minHt, col='green',
         border='white')
    overlap <- sum(minHt * diff(histinfo$breaks))
    # Add curve if prior is a function
    if (is.function(prior))
      lines(histinfo$mids, priorInfo, lwd=2, col='brown')
    # Add text
    text(mean(breaks), 0, paste0("overlap = ", round(overlap*100), "%"), pos=3, cex=cex)

    return(overlap)
  }

library(HDInterval)
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

library(rstan)
library(bayesplot)
library(ggplot2)
library(dautility)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m2 <- stan_model(file = "Regression/changepoint.stan")
data2 <- list(N=100, D=c(rnorm(50, 100, 50), rnorm(50, 200, 50)))
f2 <- sampling(m2, data=data2, iter = 500)
f2arr <- as.array(f2)
f2df <- as.data.frame(f2)
mcmc_areas(f2arr[,1,grep("^mu*", names(f2), value = T)])
mcmc_areas(f2arr[,1,grep("^sigma*", names(f2), value = T)])

saveRDS(m2, "Regression/changepoint.rds")
