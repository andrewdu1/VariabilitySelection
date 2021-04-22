# 3. cross-correlation analyses & plots

# Author: Andrew Du
# Date: 4-21-21


## Read in CMR results
CMR.res <- readRDS("CMR files/CMR results.rds")

## Read in climate variability data
clim.var <- read.csv(file = "original datasets/climate variability 250ka bins.csv", header = TRUE)

## remove time bins from climate variability that aren't in CMR results
clim.var <- clim.var[seq(1, which(clim.var$age.bin_ka == "3625")), ]

rownames(clim.var) <- clim.var$age.bin_ka
clim.var <- clim.var[, -1]

## run CCF analyses

ccf.res <- lapply(CMR.res, function(cmr){ # iterate through each taxonomic treatment
  
  apply(clim.var, 2, function(clim){ # iterate through each climate variable
    
    ### positive lags in ccf() are those where x is younger than y
    orig <- ccf(rev(1 - cmr$estimate[grep("Gamma", rownames(cmr))]), clim[-length(clim)])[0:8]
    
    extinct <- ccf(rev(1 - cmr$estimate[grep("Phi", rownames(cmr))]), clim[-1])[0:8]
    
    return(list(orig = orig, extinct = extinct))
  })
})

## create function for calculating p-values (two-tailed) from CCF results
ccf_p.value <- function(ccf.res){
  
  ccf_hat <- as.numeric(ccf.res$acf)
  ccf_se <- 1 / sqrt(ccf.res$n.used)
  
  res <- sapply(ccf_hat, function(x) ifelse(x < 0, pnorm(x, sd = ccf_se) * 2, pnorm(x, sd = ccf_se, lower.tail = FALSE) * 2))
  
  return(res)
}

## calculate P-values for CCF analyses using ccf_p.value function
ccf_p.res <- lapply(ccf.res, function(cmr){
  
  lapply(cmr, function(clim){
    
    orig <- ccf_p.value(clim$orig)
    extinct <- ccf_p.value(clim$extinct)
    
    return(list(orig = orig, extinct = extinct))
  })
})

## see how many P-values are significant before correction for multiple comparisons
raw_p.vals <- unlist(ccf_p.res) # 7 out of 180 comparisons are significant

p.res_BH <- p.adjust(raw_p.vals, method = "BH") # no comparisons are significant after BH correction



####################################################

# Plot CCF results (Fig. S9)

layout(matrix(c(1, 2, 11, 12,
                3, 4, 13, 14,
                5, 6, 15, 16,
                7, 8, 17, 18,
                9, 10, 19, 20), nrow = 5, ncol = 4, byrow = TRUE))

for(taxon in seq_along(ccf.res)){
  
  ccf.res1 <- ccf.res[[taxon]]
  
  for(clim in seq_along(ccf.res1)){
    
    clim1 <- ccf.res1[[clim]]
    
    plot(clim1$orig, ci = 0, lwd = 2, ylim = c(-1, 1), ylab = "Cross-correlation", main = paste(names(ccf.res1)[clim], "(origin.)"), cex.axis = 1.5, cex.lab = 1.5)
    
    plot(clim1$extinct, ci = 0, lwd = 2, ylim = c(-1, 1), ylab = "Cross-correlation", main = paste(names(ccf.res1)[clim], "(extinct.)"), cex.axis = 1.5, cex.lab = 1.5)
  }
}