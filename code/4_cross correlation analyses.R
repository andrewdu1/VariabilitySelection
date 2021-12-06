# 3. cross-correlation analyses & plots

# Author: Andrew Du
# Date: 4-21-21 (revised 7-27-21)


## Read in CMR results
CMR.res <- readRDS("CMR files/CMR results.rds")

turk_CMR.res <- readRDS("CMR files/Turkana CMR results.rds") # Turkana only

## Read in climate variability data
clim.var <- read.csv(file = "original datasets/climate variability 250ka bins.csv", header = TRUE, row.names = 1)

## remove time bins from climate variability that aren't in CMR results
clim.var1 <- clim.var[seq(1, which(rownames(clim.var) == "3625")), ]

clim.var1 <- clim.var1[, -ncol(clim.var1)]

d13C <- clim.var$Turkana_psol_d13C # Turkana psol data only

d13C <- d13C[rownames(clim.var) %in% seq(1125, 4125, 250)] # get out data synchronous with mammal data

## run CCF analyses

# get out origination and extinction estimates for lumped cf treatment
cmr_cf.lump <- CMR.res$neo_cf.lump # all localities cf lumped

cf.lump_orig <- rev(1 - cmr_cf.lump$estimate[grep("Gamma", rownames(cmr_cf.lump))]) # all localities origination

cf.lump_extinct <- rev(1 - cmr_cf.lump$estimate[grep("Phi", rownames(cmr_cf.lump))]) # all localities extinction

turk_cmr_lump <- turk_CMR.res$neo_cf.lump # turkana cf lumped

turk_lump_orig <- rev(1 - turk_cmr_lump$estimate[grep("Gamma", rownames(turk_cmr_lump))]) # turkana origination

turk_lump_extinct <- rev(1 - turk_cmr_lump$estimate[grep("Phi", rownames(turk_cmr_lump))]) # turkana extinction

# run CCFs comparing climate and taxonomic rates
ccf.res <- apply(clim.var1, 2, function(clim){ # iterate through each climate variable
  
  ### positive lags in ccf() are those where x is younger than y
  orig.ccf <- ccf(cf.lump_orig, clim[-length(clim)], plot = FALSE)[0:2]
  
  extinct.ccf <- ccf(cf.lump_extinct, clim[-1], plot = FALSE)[0:2]
  
  return(list(orig = orig.ccf, extinct = extinct.ccf))
})

# run CCFs for Turkana and psol data (lag up to 2)
turk.ccf <- list(
  orig = ccf(turk_lump_orig, d13C[-length(d13C)], plot = FALSE)[0:2], 
  extinct = ccf(turk_lump_extinct, d13C[-1], plot = FALSE)[0:2]
)

## create function for calculating p-values (two-tailed) from CCF results
ccf_p.value <- function(ccf.res){
  
  ccf_hat <- as.numeric(ccf.res$acf)
  ccf_se <- 1 / sqrt(ccf.res$n.used)
  
  res <- sapply(ccf_hat, function(x) ifelse(x < 0, pnorm(x, sd = ccf_se) * 2, pnorm(x, sd = ccf_se, lower.tail = FALSE) * 2))
  
  return(res)
}

## calculate P-values for CCF analyses using ccf_p.value function
ccf_p.res <- lapply(ccf.res, function(clim){
  
  orig <- ccf_p.value(clim$orig)
  extinct <- ccf_p.value(clim$extinct)
  
  return(list(orig = orig, extinct = extinct))
}) # all localities

turk_p.res <- lapply(turk.ccf, ccf_p.value) # turkana only


## see how many P-values are significant
raw_p.vals <- c(unlist(ccf_p.res), unlist(turk_p.res)) # 4 out of 36 are significant

p.res_BH <- p.adjust(raw_p.vals, method = "BH") # no comparisons are significant after BH correction


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

# New CCF plot (updated: 12-5-21)
clim.bar <- lapply(ccf.res, function(x){
  
  x1 <- cbind(x$orig$acf, x$extinct$acf)
  colnames(x1) <- c("Origination", "Extinction")
  rownames(x1) <- 0:2
  
  return(x1)
})

clim.bar$Turkana_psol_d13C <- cbind(turk.ccf$orig$acf, turk.ccf$extinct$acf)
colnames(clim.bar$Turkana_psol_d13C) <- c("Origination", "Extinction")
rownames(clim.bar$Turkana_psol_d13C) <- 0:2

clim.bar <- clim.bar[order(names(clim.bar))]


par(mfrow = c(3, 2), mar = c(5, 4, 4, 2) + 0.1)

i <- 1

barplot(clim.bar[[i]], beside = TRUE, ylab = "Cross-correlation", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5, legend.text = 0:2, args.legend = list(x = "topleft", cex = 1.5, title = "Number of lags", bty = "n"), main = names(clim.bar)[i], cex.main = 2, ylim = c(-1, 1))
abline(h = 0)
abline(v = 4.5, lty = 2)

for(i in seq_along(clim.bar)[-1]){
  
  barplot(clim.bar[[i]], beside = TRUE, ylab = "Cross-correlation", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5, main = names(clim.bar)[i], cex.main = 2, ylim = c(-1, 1))
  abline(h = 0)
  abline(v = 4.5, lty = 2)
}
