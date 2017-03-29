##' Produce a presence/absence matrix from a matrix of shared species between sites.
##' @param sharedmatrix a matrix where each cell contains the number of species shared between the column and row.
##' @param nsamples number of times to sample the possible new rows. Higher values take longer, but they are more likely to return the actual answers.
##' @return Returns a presence/absence matrix with sites across the rows and species down the columns.
##' @details This function is designed to link alpha and beta diversity values to presence/absence for maximum coverage.
##' @details Each cell in the input matrix should be an integer number of species shared between sites. So the diagonal will be alpha diversities while other cells will represent a type of beta diversity.
##' @export
produce_presabs_matrix <- function(sharedmatrix,nsamples = 1000) {
  test <- sample_newrow_presabs(sharedmatrix = sharedmatrix,nsamples = nsamples)
  num.columns<- lapply(test,ncol)
  test <-test[num.columns == max(unlist(num.columns))]
  for (i in 1:(nrow(sharedmatrix)-3)) {
    test <- sample_newrow_presabs(sharedmatrix = sharedmatrix,startingmatrix = test[[1]],nsamples = nsamples)
    if (class(test) == "list") {
      num.columns<- lapply(test,ncol)
      test <-test[num.columns == max(unlist(num.columns))] #picks the returned matrix with the most columns
    }
  }
  
  final <- test[[1]]
  colnames(final) <- paste0(1:ncol(final))
  rownames(final) <- rownames(sharedmatrix)
  final
}
