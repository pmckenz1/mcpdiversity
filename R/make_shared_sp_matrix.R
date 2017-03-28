##' Produce a "number of shared species" matrix from a presence/absence matrix
##' @param pres_abs_matrix a presence/absence matrix with sites as rows and species as columns
##' @return matrix where each cell contains the number of species shared between the column and row.
##' @details this might be useful for testing the reconstruction of a presence/absence matrix from a "shared species" matrix that includes alpha and pairwise beta values for sites.
##' @export

make_shared_sp_matrix <- function(pres_abs_matrix) {
numrows <- nrow(pres_abs_matrix)
num_shared_species_matrix <- matrix(nrow=numrows,ncol=numrows)
for (i in 1:numrows) {
  for (q in 1:numrows) {
    counting_shared <- rbind(pres_abs_matrix[i,],pres_abs_matrix[q,])
    num_shared_species_matrix[i,q] <- sum(colSums(counting_shared) > 1)
  }
}
num_shared_species_matrix
}
