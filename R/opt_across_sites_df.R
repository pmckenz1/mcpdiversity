##' Select the sites that will maximize the number of present species
##' @param numbersites the number of sites to select
##' @param data a presence/absence matrix with sites as rows and species as columns
##' @param destination directory, including file name with .csv suffix, to which a .csv of solutions will be saved
##' @return sorted list of selected sites with accompanying scores
##' @details This solves the maximum coverage problem. Given a matrix of presence/absence data for all sites, this picks the best combination of a user-defined-size subset of those sites to maximize the amount of "presence" across the selected sites.
##' @details WARNING: Because this considers every possible combination in a wildly inefficient way, it gets out-of-hand very quickly. For example, picking the best 6 sites out of 27 total requires ~296000 computations.
##' @export
opt_across_sites_df <- function(numbersites, data, destination = NULL) {
  vectors <- list()
  for (i in 1:choose(nrow(data),numbersites)) { # for 1 in the number of ways to choose # sites specified from total number of sites
    zeroes <- rep(0,nrow(data))                 # vector of zeroes length of the total number of sites
    zeroes[combn(1:nrow(data),numbersites)[,i]] <- 1 # inserts 1 for each way to pick the # sites specified
    vectors[[i]] <- as.logical(zeroes)  #builds list of ways to choose # sites specified
    print(paste0(i,"/",choose(nrow(data),numbersites)," combinations calculated."))
  }

  notsorted_mainlist <- list()
  for (i in 1:length(vectors)) {
    temp.selected <- data[vectors[[i]],]
    num.species <- sum(colSums(temp.selected) > 0) # total number species present
    percent.present <- num.species/ncol(data) # percent of total species present across selected sites
    selected.sites <- row.names(temp.selected) # vector of selected site names
    species.present <- colnames(data)[colSums(temp.selected) > 0] # vector of species names
    templist <- list("num.species" = num.species, "percent.present" = percent.present,
                     "selected.sites" = selected.sites, "species.present" = species.present)
    notsorted_mainlist[[i]] <- templist
    print(paste0(i,"/",length(vectors)," combinations summed."))
  }
  print("Sorting...")
  mainlist <- notsorted_mainlist[order(unlist(lapply(notsorted_mainlist,function(x) head(x, 1))),decreasing = TRUE)]
  #Now write out the results as a data frame
  df <- t(data.frame(lapply(mainlist,function(x) x[3])))
  df2 <- t(data.frame(lapply(mainlist,function(x) x[c(1)])))
  df3 <- t(data.frame(lapply(mainlist,function(x) x[c(2)])))
  fulldf <- cbind.data.frame(df,df2,df3)
  colheads<-sapply("species",paste0,1:numbersites)
  colnames(fulldf) <- c(colheads,"num.species","percent.present")
  if (!is.null(destination)) {
    write.csv(fulldf,destination,row.names = TRUE)
  }
  mainlist
}

