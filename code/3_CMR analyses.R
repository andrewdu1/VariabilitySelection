# 3. run & plot capture-mark-recapture (CMR) analyses

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

## save results
saveRDS(CMR.res, file = "CMR files/CMR results.rds")



####################################################

# Plotting the CMR results. NB: plots were modified in Illustrator after being created in R

## create time bins of 0.25 Myr from 3.75-0 Ma
bins <- seq(3.75, 0, -0.25)

## bin age midpoints for plotting
time.bin.plot <- rev(bins)[-length(bins)] + .25 / 2

## (1) cf. lumped results
CMR_cf.lump <- CMR.res$neo_cf.lump

### Origination probability
plot(time.bin.plot[-length(time.bin.plot)], 
     rev(1 - CMR_cf.lump$estimate[grep("Gamma", rownames(CMR_cf.lump))]), 
     ylim = c(0, 1), 
     type = "l",
     lwd = 3,
     col = "blue", 
     xlim = range(bins), 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot[-length(time.bin.plot)], rev(time.bin.plot[-length(time.bin.plot)])),
        c(rev(1 - CMR_cf.lump$lcl[grep("Gamma", rownames(CMR_cf.lump))]), 1 - CMR_cf.lump$ucl[grep("Gamma", rownames(CMR_cf.lump))]),
        col = adjustcolor("blue", alpha.f = 0.2),
        border = NA)

### Add extinction probability
lines(time.bin.plot[-1], 
      rev(1 - CMR_cf.lump$estimate[grep("Phi", rownames(CMR_cf.lump))]), 
      lwd = 3,
      col = "red")
polygon(c(time.bin.plot[-1], rev(time.bin.plot[-1])),
        c(rev(1 - CMR_cf.lump$lcl[grep("Phi", rownames(CMR_cf.lump))]), 1 - CMR_cf.lump$ucl[grep("Phi", rownames(CMR_cf.lump))]),
        col = adjustcolor("red", alpha.f = 0.2),
        border = NA)

### Sampling probability
plot(time.bin.plot, 
     rev(CMR_cf.lump$estimate[grep("p", rownames(CMR_cf.lump))]), 
     ylim = c(0, 1), 
     type = "l",
     lty = 2, 
     lwd = 3,
     xlim = range(bins), 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot, rev(time.bin.plot)),
        c(rev(CMR_cf.lump$lcl[grep("p", rownames(CMR_cf.lump))]), CMR_cf.lump$ucl[grep("p", rownames(CMR_cf.lump))]),
        col = adjustcolor("gray", alpha.f = 0.2),
        border = NA)


## (2) dropping cf results
CMR_drop <- CMR.res$neo_drop

### Origination probability
plot(time.bin.plot[-length(time.bin.plot)], 
     rev(1 - CMR_drop$estimate[grep("Gamma", rownames(CMR_drop))]), 
     ylim = c(0, 1), 
     type = "l",
     lwd = 3,
     col = "blue", 
     xlim = range(bins), 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot[-length(time.bin.plot)], rev(time.bin.plot[-length(time.bin.plot)])),
        c(rev(1 - CMR_drop$lcl[grep("Gamma", rownames(CMR_drop))]), 1 - CMR_drop$ucl[grep("Gamma", rownames(CMR_drop))]),
        col = adjustcolor("blue", alpha.f = 0.2),
        border = NA)

### Add extinction probability
lines(time.bin.plot[-1], 
      rev(1 - CMR_drop$estimate[grep("Phi", rownames(CMR_drop))]), 
      lwd = 3,
      col = "red")
polygon(c(time.bin.plot[-1], rev(time.bin.plot[-1])),
        c(rev(1 - CMR_drop$lcl[grep("Phi", rownames(CMR_drop))]), 1 - CMR_drop$ucl[grep("Phi", rownames(CMR_drop))]),
        col = adjustcolor("red", alpha.f = 0.2),
        border = NA)

### Sampling probability
plot(time.bin.plot, 
     rev(CMR_drop$estimate[grep("p", rownames(CMR_drop))]), 
     ylim = c(0, 1), 
     type = "l",
     lty = 2, 
     lwd = 3,
     xlim = range(bins), 
     xlab = "Time (Ma)", 
     ylab = "Probability", 
     cex.axis = 1.5, 
     cex.lab = 1.5)
polygon(c(time.bin.plot, rev(time.bin.plot)),
        c(rev(CMR_drop$lcl[grep("p", rownames(CMR_drop))]), CMR_drop$ucl[grep("p", rownames(CMR_drop))]),
        col = adjustcolor("gray", alpha.f = 0.2),
        border = NA)