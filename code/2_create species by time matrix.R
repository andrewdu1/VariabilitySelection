# 2. create species by time matrix

# Author: Andrew Du
# Date: 4-21-21 (revised 7-27-21)


## read in PA data
pa_cf.lump <- read.csv("modified datasets/PA matrix lump cf.csv", header = TRUE, row.names = 1)
pa_drop <- read.csv("modified datasets/PA matrix drop cf.csv", header = TRUE, row.names = 1)

## read in site age data
site <- read.csv("original datasets/site data.csv", header = TRUE, strip.white = TRUE)

## get out Turkana PA data only
turk.site <- site[site$Region == "Turkana", ]

turk_pa_cf.lump <- pa_cf.lump[, colnames(pa_cf.lump) %in% turk.site$Unit]
turk_pa_drop <- pa_drop[, colnames(pa_drop) %in% turk.sites$Unit]

## create time bins of 0.25 Myr from 3.75-0 Ma. For Turkana, 4.25-1 Ma
bins <- seq(3.75, 0, -0.25)
turk.bins <- seq(4.25, 1, -0.25)

site1 <- data.frame(site, time.bin = cut(site$Mean.Age, bins))
turk.site <- data.frame(turk.site, time.bin = cut(turk.site$Mean.Age, turk.bins))

## create species by time matrices
time.bins <- sort(unique(site1$time.bin), decreasing = TRUE) # vector of time bins

spec.time_cf.lump <- matrix(data = NA, nrow = nrow(pa_cf.lump), ncol = length(time.bins), dimnames = list(rownames(pa_cf.lump), time.bins))

spec.time_drop <- matrix(data = NA, nrow = nrow(pa_drop), ncol = length(time.bins), dimnames = list(rownames(pa_drop), time.bins))

for(i in seq_along(time.bins)){
  
  time.bin <- time.bins[i] # work with one time bin at a time
  
  sites.in.bin <- site1$Unit[site1$time.bin == time.bin] # get out sites belong to time bin
  
  pa.bin_cf.lump <- pa_cf.lump[, colnames(pa_cf.lump) %in% sites.in.bin] # subset PA matrix for that time bin
  
  pa.bin_drop <- pa_drop[, colnames(pa_drop) %in% sites.in.bin] # subset PA matrix for that time bin
  
  if(is.null(dim(pa.bin_cf.lump))){ # in case there's only one site in the time bin
    
    spec.time_cf.lump[, i] <- pa.bin_cf.lump
    spec.time_drop[, i] <- pa.bin_drop
    
  } else {
    
    spec.time_cf.lump[, i] <- rowSums(pa.bin_cf.lump)
    spec.time_drop[, i] <- rowSums(pa.bin_drop)
  }
}

## change all numbers greater than zero to one
spec.time_cf.lump[spec.time_cf.lump > 0] <- 1
spec.time_drop[spec.time_drop > 0] <- 1

## remove species not found in any time bin
spec.time_cf.lump <- spec.time_cf.lump[rowSums(spec.time_cf.lump) > 0, ]
spec.time_drop <- spec.time_drop[rowSums(spec.time_drop) > 0, ]


## create species by time matrices for Turkana
turk.time.bins <- sort(unique(turk.site$time.bin), decreasing = TRUE) # vector of time bins

turk_spec.time_cf.lump <- matrix(data = NA, nrow = nrow(turk_pa_cf.lump), ncol = length(turk.time.bins), dimnames = list(rownames(turk_pa_cf.lump), turk.time.bins))

turk_spec.time_drop <- matrix(data = NA, nrow = nrow(turk_pa_drop), ncol = length(turk.time.bins), dimnames = list(rownames(turk_pa_drop), turk.time.bins))

for(i in seq_along(turk.time.bins)){
  
  time.bin <- turk.time.bins[i] # work with one time bin at a time
  
  sites.in.bin <- turk.site$Unit[turk.site$time.bin == time.bin] # get out sites belong to time bin
  
  turk_pa.bin_cf.lump <- turk_pa_cf.lump[, colnames(turk_pa_cf.lump) %in% sites.in.bin] # subset PA matrix for that time bin
  
  turk_pa.bin_drop <- turk_pa_drop[, colnames(turk_pa_drop) %in% sites.in.bin] # subset PA matrix for that time bin
  
  if(is.null(dim(turk_pa.bin_cf.lump))){ # in case there's only one site in the time bin
    
    turk_spec.time_cf.lump[, i] <- turk_pa.bin_cf.lump
    turk_spec.time_drop[, i] <- turk_pa.bin_drop
    
  } else {
    
    turk_spec.time_cf.lump[, i] <- rowSums(turk_pa.bin_cf.lump)
    turk_spec.time_drop[, i] <- rowSums(turk_pa.bin_drop)
  }
}

## change all numbers greater than zero to one
turk_spec.time_cf.lump[turk_spec.time_cf.lump > 0] <- 1
turk_spec.time_drop[turk_spec.time_drop > 0] <- 1

## remove species not found in any time bin
turk_spec.time_cf.lump <- turk_spec.time_cf.lump[rowSums(turk_spec.time_cf.lump) > 0, ]
turk_spec.time_drop <- turk_spec.time_drop[rowSums(turk_spec.time_drop) > 0, ]


## write out species-by-time matrices
write.csv(spec.time_cf.lump, "modified datasets/species by time matrix lump cf.csv")
write.csv(spec.time_drop, "modified datasets/species by time matrix drop cf.csv")

write.csv(turk_spec.time_cf.lump, "modified Datasets/Turkana species by time matrix lump cf.csv")
write.csv(turk_spec.time_drop, "modified Datasets/Turkana species by time matrix drop cf.csv")