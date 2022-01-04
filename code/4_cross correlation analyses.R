# 4. cross-correlation analyses & plots

# Author: Andrew Du
# Date: 4-21-21 (revised 1-3-21)


## Read in CMR results
CMR.res <- readRDS("CMR files/CMR results.rds")

turk_CMR.res <- readRDS("CMR files/Turkana CMR results.rds") # Turkana only

# get out origination and extinction estimates for lumped cf treatment & combine into one list
cmr_cf.lump <- CMR.res$neo_cf.lump # all localities cf lumped
turk_cmr_lump <- turk_CMR.res$neo_cf.lump # turkana cf lumped

cmr.rates <- list(
  cf.lump_orig = 1 - cmr_cf.lump$estimate[grep("Gamma", rownames(cmr_cf.lump))], # all localities origination
  cf.lump_extinct = 1 - cmr_cf.lump$estimate[grep("Phi", rownames(cmr_cf.lump))], # all localities extinction
  turk_lump_orig = 1 - turk_cmr_lump$estimate[grep("Gamma", rownames(turk_cmr_lump))], # turkana origination
  turk_lump_extinct = 1 - turk_cmr_lump$estimate[grep("Phi", rownames(turk_cmr_lump))] # turkana extinction
)

# create a list of ages, corresponding to rates in cmr.rates
bins <- seq(3.75, 0, -0.25) # create time bins of 0.25 Myr from 3.75-0 Ma
time.bin.mid <- bins[-length(bins)] + diff(bins)[1] / 2 # bin age midpoints

turk.bin.mid <- seq(4.125, 1.125, -0.25) # time bin midpoints for 0.25 Myr from 4.125-1.125 for Turkana

ages <- list(
  cf.lump_orig_ages = time.bin.mid[-1],
  cf.lump_extinct_ages = time.bin.mid[-length(time.bin.mid)],
  turk_lump_orig_ages = turk.bin.mid[-1],
  turk_lump_extinct_ages = turk.bin.mid[-length(turk.bin.mid)]
)

## Read in climate variability data
clim.var <- read.csv(file = "original datasets/climate variability 250ka bins.csv", header = TRUE, row.names = 1)

## remove time bins from climate variability that aren't in CMR results
clim.var1 <- clim.var[seq(1, which(rownames(clim.var) == "3625")), ]

clim.var1 <- clim.var1[, -ncol(clim.var1)]

clim.var1 <- clim.var1[seq(nrow(clim.var1), 1), ] # reverse row order, so time bins match those of CMR rates

d13C <- clim.var$Turkana_psol_d13C # Turkana psol data only

d13C <- d13C[rownames(clim.var) %in% seq(1125, 4125, 250)] # get out data synchronous with mammal data

d13C <- rev(d13C) # reverse row order, so time bins match those of CMR rates


####################################################

## Detrend each time series using a LOWESS regression
# write a function to do this
  # ARGUMENTS:
    # age: age associated with each data point in the time series
    # x: variable value associated with each data point
    # span: smoothing span of LOWESS regression (passed on to lowess() function)
detrend <- function(age, x, span = 2 / 3){
  
  lowess.res <- lowess(age, x, f = span)
  return(x - lowess.res$y)
}

cmr.rates.detrend <- mapply(detrend, age = ages, x = cmr.rates)


## run CCF analyses
# create function for estimate CCF and its p-value
  # ARGUMENTS: 
    # climate: climate time series
    # cmr_rate: CMR rate time series
    # ccf.lag: max #lags to consider in CCF 
    # acf.lag: max. # lags in ACF for computing SE 
CCF <- function(climate, cmr_rate, ccf.lag, acf.lag){
  
  n <- length(climate)
  
  ccf.res <- ccf(cmr_rate, climate, plot = FALSE)[0:ccf.lag]
  
  ccf.hat <- ccf.res$acf
  
  acf.clim <- acf(climate, lag.max = acf.lag, plot = FALSE)$acf
  acf.cmr <- acf(cmr_rate, lag.max = acf.lag, plot = FALSE)$acf
  
  ccf.se <- sqrt((1 + 2 * sum(acf.clim * acf.cmr)) / n)
  
  p.vals <- sapply(ccf.hat, function(x) pnorm(abs(x), sd = ccf.se, lower.tail = FALSE) * 2)
  
  return(list(lag.n = seq(0, ccf.lag), ccf.hat = ccf.hat, ccf.se = ccf.se, p.vals = p.vals))
}

# run CCFs for E. African data
ccf.res.orig <- apply(clim.var1, 2, function(clim) CCF(climate = clim[-1], cmr_rate = cmr.rates$cf.lump_orig, ccf.lag = 2, acf.lag = length(clim) - 1))

ccf.res.extinct <- apply(clim.var1, 2, function(clim) CCF(climate = clim[-length(clim)], cmr_rate = cmr.rates$cf.lump_extinct, ccf.lag = 2, acf.lag = length(clim) - 1))

# run CCFs for Turkana data
ccf.turk.orig <- CCF(climate = d13C[-1], cmr_rate = cmr.rates$turk_lump_orig, ccf.lag = 2, acf.lag = length(d13C) - 1)

ccf.turk.extinct <- CCF(climate = d13C[-length(d13C)], cmr_rate = cmr.rates$turk_lump_extinct, ccf.lag = 2, acf.lag = length(d13C) - 1)

## see how many P-values are significant
p.vals <- c(sapply(ccf.res.orig, function(x) x$p.vals), 
            sapply(ccf.res.extinct, function(x) x$p.vals),
            ccf.turk.orig$p.vals,
            ccf.turk.extinct$p.vals)

sum(p.vals <= 0.05) # none are significant


####################################################

# create table of results
ccf.hat <- sapply(ccf.res, function(clim) sapply(clim, function(x) x$acf))

turk.hat <- sapply(turk.ccf, function(x) x$acf)

ccf.df <- data.frame(ccf = c(c(ccf.hat), c(turk.hat)), 
                        clim.var = rep(colnames(clim.var), each = 6), 
                        rate = rep(rep(c("origination", "extinction"), each = 3), 6), 
                        lag = rep(0:2, 12),
                        raw.p = raw_p.vals,
                        BH.p = p.res_BH)


####################################################

# CCF plot
par(mfrow = c(3, 2), mar = c(5 - 1, 4, 4, 2) + 0.1)

# plot the first CCF results by itself, so I can add a legend
i <- 1

ccf.mat <- cbind(as.numeric(ccf.res.orig[[i]]$ccf.hat), as.numeric(ccf.res.extinct[[i]]$ccf.hat))

barplot(ccf.mat, beside = TRUE, ylim = c(-1, 1), ylab = "Cross-correlation", names.arg = c("Origination", "Extinction"), main = names(ccf.res.orig)[i], cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5, cex.main = 1.5, legend.text = ccf.res.orig[[i]]$lag.n, args.legend = list(x = 3.5, y = 1, title = "Number of lags", bty = "n", cex = 1.25))

abline(h = 0)
abline(v = 4.5)

# origination significance thresholds
segments(x0 = 0, x1 = 4.5, y0 = ccf.res.orig[[i]]$ccf.se * qnorm(0.025), lty = 2)
segments(x0 = 0, x1 = 4.5, y0 = ccf.res.orig[[i]]$ccf.se * qnorm(0.975), lty = 2)

# extinction significance thresholds
segments(x0 = 4.5, x1 = 9, y0 = ccf.res.extinct[[i]]$ccf.se * qnorm(0.025), lty = 2)
segments(x0 = 4.5, x1 = 9, y0 = ccf.res.extinct[[i]]$ccf.se * qnorm(0.975), lty = 2)

# iterate through rest of CCF results
for(i in seq_along(ccf.res.orig)[-1]){
  
  ccf.mat <- cbind(as.numeric(ccf.res.orig[[i]]$ccf.hat), as.numeric(ccf.res.extinct[[i]]$ccf.hat))
  
  barplot(ccf.mat, beside = TRUE, ylim = c(-1, 1), ylab = "Cross-correlation", names.arg = c("Origination", "Extinction"), main = names(ccf.res.orig)[i], cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5, cex.main = 1.5)
  
  abline(h = 0)
  abline(v = 4.5)
  
  # origination significance thresholds
  segments(x0 = 0, x1 = 4.5, y0 = ccf.res.orig[[i]]$ccf.se * qnorm(0.025), lty = 2)
  segments(x0 = 0, x1 = 4.5, y0 = ccf.res.orig[[i]]$ccf.se * qnorm(0.975), lty = 2)
  
  # extinction significance thresholds
  segments(x0 = 4.5, x1 = 9, y0 = ccf.res.extinct[[i]]$ccf.se * qnorm(0.025), lty = 2)
  segments(x0 = 4.5, x1 = 9, y0 = ccf.res.extinct[[i]]$ccf.se * qnorm(0.975), lty = 2)
}

# Turkana CCF results
barplot(cbind(as.numeric(ccf.turk.orig$ccf.hat), as.numeric(ccf.turk.extinct$ccf.hat)), beside = TRUE, ylim = c(-1, 1), ylab = "Cross-correlation", main = colnames(clim.var)[ncol(clim.var)], names.arg = c("Origination", "Extinction"), cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5, cex.main = 1.5)

abline(h = 0)
abline(v = 4.5)

# origination significance thresholds
segments(x0 = 0, x1 = 4.5, y0 = ccf.turk.orig$ccf.se * qnorm(0.025), lty = 2)
segments(x0 = 0, x1 = 4.5, y0 = ccf.turk.extinct$ccf.se * qnorm(0.975), lty = 2)

# extinction significance thresholds
segments(x0 = 4.5, x1 = 9, y0 = ccf.res.extinct[[i]]$ccf.se * qnorm(0.025), lty = 2)
segments(x0 = 4.5, x1 = 9, y0 = ccf.res.extinct[[i]]$ccf.se * qnorm(0.975), lty = 2)