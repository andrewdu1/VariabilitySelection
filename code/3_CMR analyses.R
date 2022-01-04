# 3. run & plot capture-mark-recapture (CMR) analyses

# Author: Andrew Du
# Date: 4-21-21 (revised 1-3-22)


## read in species-by-time matrices
spec.time_cf.lump <- read.csv("modified datasets/species by time matrix lump cf.csv", header = TRUE, row.names = 1)
spec.time_drop <- read.csv("modified datasets/species by time matrix drop cf.csv", header = TRUE, row.names = 1)

turk_spec.time_cf.lump <- read.csv("modified Datasets/Turkana species by time matrix lump cf.csv", header = TRUE, row.names = 1)
turk_spec.time_drop <- read.csv("modified Datasets/Turkana species by time matrix drop cf.csv", header = TRUE, row.names = 1)

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

## save results
saveRDS(CMR.res, file = "CMR files/CMR results.rds")


### Do the same but for Turkana
## Convert matrices into .inp files
turk_filename_cf.lump <- "CMR files/Turkana_dataMatrix_CMR_cfLump.inp"
turk_filename_drop <- "CMR files/Turkana_dataMatrix_CMR_drop.inp"

## Code modified after Smiley 2018 in Paleobiology
for(i in seq_len(nrow(turk_spec.time_cf.lump))){
  vec <- paste(c(turk_spec.time_cf.lump[i, ], " 1;"), collapse = "")
  write(vec, turk_filename_cf.lump, append = ifelse(i == 1, FALSE, TRUE))
}

for(i in seq_len(nrow(turk_spec.time_drop))){
  vec <- paste(c(turk_spec.time_drop[i, ], " 1;"), collapse = "")
  write(vec, turk_filename_drop, append = ifelse(i == 1, FALSE, TRUE))
}

turk_neo_cf.lump <- convert.inp(turk_filename_cf.lump)
turk_neo_drop <- convert.inp(turk_filename_drop)

turk_neo <- list(neo_cf.lump = turk_neo_cf.lump, neo_drop = turk_neo_drop)

## process data to fit Pradsen model
turk_neo.prad <- lapply(turk_neo, function(x) process.data(x, model = "Pradsen"))

## run time-varying model
turk_CMR <- lapply(turk_neo.prad, function(x) mark(x, model.parameters = list(Phi = Phi.t, p = p.t, Gamma = Gamma.t), delete = TRUE))

## subset out estimated parameters of interest
turk_CMR.res <- lapply(turk_CMR, function(x) x$results$real)

## save results
saveRDS(turk_CMR.res, file = "CMR files/Turkana CMR results.rds")


####################################################

# Plotting the CMR results (did not do this for Turkana results, but the code is the same). NB: plots were modified in Illustrator after being created in R

## create time bins of 0.25 Myr from 3.75-0 Ma
bins <- seq(3.75, 0, -0.25)

## bin age midpoints for plotting
time.bin.plot <- bins[-length(bins)] - .25 / 2

## (1) cf. lumped results
CMR_cf.lump <- CMR.res$neo_cf.lump

### Origination probability (1-gamma_i: time i-1 to i)
plot(time.bin.plot[-1], 
     1 - CMR_cf.lump$estimate[grep("Gamma", rownames(CMR_cf.lump))], 
     xlim = range(bins),
     ylim = c(0, 1), 
     type = "l",
     lwd = 3,
     col = "blue", 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot[-1], rev(time.bin.plot[-1])),
        c(1 - CMR_cf.lump$lcl[grep("Gamma", rownames(CMR_cf.lump))], rev(1 - CMR_cf.lump$ucl[grep("Gamma", rownames(CMR_cf.lump))])),
        col = adjustcolor("blue", alpha.f = 0.2),
        border = NA)

### Add extinction probability (1-phi_i: time i to i+1)
lines(time.bin.plot[-length(time.bin.plot)], 
      1 - CMR_cf.lump$estimate[grep("Phi", rownames(CMR_cf.lump))], 
      lwd = 3,
      col = "red")
polygon(c(time.bin.plot[-length(time.bin.plot)], rev(time.bin.plot[-length(time.bin.plot)])),
        c(1 - CMR_cf.lump$lcl[grep("Phi", rownames(CMR_cf.lump))], rev(1 - CMR_cf.lump$ucl[grep("Phi", rownames(CMR_cf.lump))])),
        col = adjustcolor("red", alpha.f = 0.2),
        border = NA)

### Sampling probability (p_i)
plot(time.bin.plot, 
     CMR_cf.lump$estimate[grep("p", rownames(CMR_cf.lump))], 
     xlim = range(bins),
     ylim = c(0, 1), 
     type = "l",
     lty = 2, 
     lwd = 3, 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot, rev(time.bin.plot)),
        c(CMR_cf.lump$lcl[grep("p", rownames(CMR_cf.lump))], rev(CMR_cf.lump$ucl[grep("p", rownames(CMR_cf.lump))])),
        col = adjustcolor("gray", alpha.f = 0.2),
        border = NA)


## (2) dropping cf results
CMR_drop <- CMR.res$neo_drop

### Origination probability
plot(time.bin.plot[-1], 
     1 - CMR_drop$estimate[grep("Gamma", rownames(CMR_drop))], 
     xlim = range(bins),
     ylim = c(0, 1), 
     type = "l",
     lwd = 3,
     col = "blue", 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot[-1], rev(time.bin.plot[-1])),
        c(1 - CMR_drop$lcl[grep("Gamma", rownames(CMR_drop))], rev(1 - CMR_drop$ucl[grep("Gamma", rownames(CMR_drop))])),
        col = adjustcolor("blue", alpha.f = 0.2),
        border = NA)

### Add extinction probability
lines(time.bin.plot[-length(time.bin.plot)], 
      1 - CMR_drop$estimate[grep("Phi", rownames(CMR_drop))], 
      lwd = 3,
      col = "red")
polygon(c(time.bin.plot[-length(time.bin.plot)], rev(time.bin.plot[-length(time.bin.plot)])),
        c(1 - CMR_drop$lcl[grep("Phi", rownames(CMR_drop))], rev(1 - CMR_drop$ucl[grep("Phi", rownames(CMR_drop))])),
        col = adjustcolor("red", alpha.f = 0.2),
        border = NA)

### Sampling probability
plot(time.bin.plot, 
     CMR_drop$estimate[grep("p", rownames(CMR_drop))], 
     xlim = range(bins),
     ylim = c(0, 1), 
     type = "l",
     lty = 2, 
     lwd = 3,
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot, rev(time.bin.plot)),
        c(CMR_drop$lcl[grep("p", rownames(CMR_drop))], rev(CMR_drop$ucl[grep("p", rownames(CMR_drop))])),
        col = adjustcolor("gray", alpha.f = 0.2),
        border = NA)