##' Produce rarity-weighted richness values for each site.
##' @param data a presence/absence matrix with sites across the rows and species down the columns. Should include site names for the rows.
##' @param nselect number of sites to select to maximize the rarity-weighted richness score.
##' @return Returns a list including a vector of the selected sites, an integer number of species included, and a vector of rarity-weighted richness scores for the selected sites.
##' @details This function selects sites based on rarity-weighted richness.
##' @export
compute_rwr <- function(data,nselect) {
### RWR - rarity weighted richness ala Albuquerque and Beier 2016 ###
  # code based on that in a script by Heather Jackson
# first calculate number of sites in which each species occurs (i.e. sum across rows)
data <- as.matrix(data)
if (is.null(rownames(data))) {
  stop("Your row names need to be labeled!")
}
nsites <- apply(data,2,sum)

# remove species (rows) that are not present in any site
data2 <- data[,!(nsites == 0)]
nsites2 <- nsites[!(nsites == 0)]

# calculate rarity (1/nsites)
rarity <- 1/nsites2

# next calculate rarity weighted richness (sum rarity scores of species within a site)
raritymat <- data2 %*% diag(rarity) #this multiplies columns by the rarity of that species
rwr <- apply(raritymat,1,sum)
rwrdf <- cbind.data.frame(rownames(data),rwr)
rwrdf <- rwrdf[order(as.numeric(as.character(rwrdf[,2])),decreasing = TRUE),]
rwrdf_selected <- rwrdf[1:nselect,]

num.species <- sum(colSums(data[as.character(rwrdf_selected$`rownames(data)`),])>0)

list(site.names = as.character(rwrdf_selected$`rownames(data)`), num.species = num.species,rwr.scores = rwrdf_selected$rwr)
}