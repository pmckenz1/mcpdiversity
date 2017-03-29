##' Produce a random presence/absence matrix of user-specified size
##' @param nrow the number of sites to select
##' @param ncol a presence/absence matrix with sites as rows and species as columns
##' @return Returns a random presence/absence matrix of the specified size.
##' @export
simulate_presabs_matrix <- function(nrow,ncol) {
sim_matrix <- matrix(nrow = nrow,ncol = ncol)
row.names(sim_matrix) <- 1:nrow(sim_matrix)
for (i in 1:nrow(sim_matrix)) {
  pres_sample <- sample((ncol(sim_matrix)/4),1)
  cols <- sample(ncol(sim_matrix),pres_sample,replace = FALSE)
  sim_matrix[i,cols] <- 1
  sim_matrix[i,!(1:ncol(sim_matrix) %in% cols)] <- 0
}
sim_matrix
}