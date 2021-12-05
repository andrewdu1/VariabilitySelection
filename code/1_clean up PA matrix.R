# 1. clean up PA matrix

# Author: Andrew Du
# Date: 4-21-21 (revised 12/5/21)


# Cleaning is going to be done in two ways:
# (1) dropping sp. and indet. species and lumping cf. species
# (2) dropping sp., indet., and cf. species

## Read in data
pa <- read.csv("original datasets/PA matrix.csv", header = TRUE, strip.white = TRUE)

# (1) lump cf species

## get rid of sp & indet
pa <- pa[pa$Qual != "sp" & pa$Qual != "indet", ]

## get rid of cf modifier
pa$Species <- as.character(pa$Species)

which.cf <- which(pa$Qual == "cf")

library(stringr)

for(i in seq_along(which.cf)) pa$Species[which.cf[i]] <- str_replace(pa$Species[which.cf[i]], "cf_", "") # replace cf with empty space

## create new PA matrix where duplicated species names are combined
pa1 <- pa[, seq(grep("Am.Ado", colnames(pa)), grep("Usno", colnames(pa)))] # get out presence-absence data only

pa.comb <- by(pa1, pa$Species, colSums) # aggregate PA data by species name

pa.comb1 <- t(simplify2array(pa.comb)) # transform to array

## remove remaining sp
which.sp <- substr(rownames(pa.comb1), nchar(rownames(pa.comb1)) - 1, nchar(rownames(pa.comb1))) == "sp"

pa_cf.lump <- pa.comb1[!which.sp, ]

# (2) drop sp., indet., and cf. species
pa2 <- pa[pa$Qual == "aff" | pa$Qual == "GOOD", ]

pa_drop <- pa2[, seq(grep("Am.Ado", colnames(pa2)), grep("Usno", colnames(pa2)))] # get out PA data only

rownames(pa_drop) <- pa2$Species

## write out cleaned PA matrices
write.csv(pa_cf.lump, "modified datasets/PA matrix lump cf.csv")
write.csv(pa_drop, "modified datasets/PA matrix drop cf.csv")