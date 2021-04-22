# 3. run capture-mark-recapture (CMR) analyses

# Author: Andrew Du
# Date: 4-21-21


## read in species-by-time matrices
spec.time_cf.lump <- read.csv("modified datasets/species by time matrix lump cf.csv", header = TRUE, row.names = 1)
spec.time_drop <- read.csv("modified datasets/species by time matrix drop cf.csv", header = TRUE, row.names = 1)

## load RMark package. NB: Need to have MARK installed on your computer first (http://www.phidot.org/software/mark/)
library(RMark)

## Convert matrices into .inp files
filename_cf.lump <- "CMR files/dataMatrix_CMR_cfLump.inp"
filename_drop <- "CMR files/dataMatrix_CMR_drop.inp"

## Code modified after Smiley 2018 in Paleobiology
for(i in seq_len(nrow(spec.time_cf.lump))){
  vec <- paste(c(spec.time_cf.lump[i, ], " 1;"), collapse = "")
  write(vec, filename_cf.lump, append = ifelse(i == 1, FALSE, TRUE))
}

for(i in seq_len(nrow(spec.time_drop))){
  vec <- paste(c(spec.time_drop[i, ], " 1;"), collapse = "")
  write(vec, filename_drop, append = ifelse(i == 1, FALSE, TRUE))
}

neo_cf.lump <- convert.inp(filename_cf.lump)
neo_drop <- convert.inp(filename_drop)

neo <- list(neo_cf.lump = neo_cf.lump, neo_drop = neo_drop)

## process data to fit Pradsen model
neo.prad <- lapply(neo, function(x) process.data(x, model = "Pradsen"))

## Model specifications (time-varying)
Phi.t <- list(formula = ~time)
Gamma.t <- list(formula = ~time)
p.t <- list(formula = ~time)

## run time-varying model
CMR <- lapply(neo.prad, function(x) mark(x, model.parameters = list(Phi = Phi.t, p = p.t, Gamma = Gamma.t), delete = TRUE))

## subset out estimated parameters of interest
CMR.res <- lapply(CMR, function(x) x$results$real)