##' Add a new row to a presence/absence matrix being reconstructed from a matrix of pairwise shared species
##' @param sharedmatrix a matrix where each cell contains the number of species shared between the column and row.
##' @param nsamples number of times to sample the possible new rows. Higher values take longer, but they are more likely to return the actual answers.
##' @param startingmatrix the first few rows (two or more) of a presence/absence matrix being reconstructed. This might be, for example, the returned value from a previous run of this function.
##' @return matrix (or matrices) of the input `startingmatrix` with an additional row of presence/absence data corresponding to the best scoring sample.
##' @details This function is the opposite of `make_shared_sp_matrix()`, which converts a presence/absence table to a matrix of pairwise shared species.
##' @details `sample_newrow_presabs()` uses a matrix of pairwise shared values to add an additional row to a presence/absence matrix being constructed from those values. If no startingmatrix is provided, this function outputs the first three rows. Otherwise, this function adds a row to the input startingmatrix.
##' @details This function randomly samples possible new rows and scores each sample based on its distance from matching the matrix of pairwise shared species.
##' @export

sample_newrow_presabs <- function(sharedmatrix,nsamples=10000,startingmatrix=NULL) {
  if (is.null(startingmatrix)) {
    shared11 <- sharedmatrix[1,1]
    shared22 <- sharedmatrix[2,2]
    shared12 <- sharedmatrix[1,2]
    nstartingcolumns <- shared11+shared22-shared12
    startingmatrix <- rbind(rep(0,nstartingcolumns),rep(0,nstartingcolumns))
    startingmatrix[1,1:shared11] <- 1
    if (shared12>0) {
      startingmatrix[2,1:shared12] <- 1
    }
    if (shared22 > shared12) {
      startingmatrix[2,(shared11+1):((shared11)+(shared22-shared12))]<- 1
    }
    currentmatrix <- startingmatrix
  }
  else {
    currentmatrix <- startingmatrix
  }
whichrow <- nrow(currentmatrix) + 1
currentmatrix <- t(apply(currentmatrix,1,c,rep(0,sharedmatrix[whichrow,whichrow])))

vector.to.row <- function(avector) {
  if (is.null(nrow(avector)) && is.vector(avector)) {
    dim(avector) <- c(1,length(avector))
  }
  avector
}


score.new.row <- function(newrow,matrix) {
  scores<-integer(0)
  for (i in 1:nrow(matrix)) {
    scores<-c(scores,abs(sharedmatrix[i,whichrow] - sum(colSums(rbind(newrow,matrix[i,]))>1)))
  }
  scores
}

scores_rows <- integer(0)
columns <- ncol(currentmatrix)
for (i in 1:nsamples) {
firstsite<-rep(0,sharedmatrix[1,1])
firstsite[sample(sharedmatrix[1,1],sharedmatrix[1,whichrow],replace=F)]<-1
othersites <- rep(0,(columns-sharedmatrix[1,1]))
othersites[sample((columns-sharedmatrix[1,1]),(sharedmatrix[whichrow,whichrow]-sharedmatrix[1,whichrow]))] <- 1

poss.new.row<-c(firstsite,othersites)
score.and.row <- c(sum(score.new.row(poss.new.row,currentmatrix)),poss.new.row)

scores_rows <- rbind(scores_rows,score.and.row)
print(paste0(i,", ncol is: ",columns))
}
#then reduce to lowest scores
scores_rows <- scores_rows[scores_rows[,1]==min(scores_rows[,1]),]

scores_rows <- vector.to.row(scores_rows) #this makes sure our object is of rows, even if there's just one...

print(paste0("Lowest row score is: ",scores_rows[1,1]))
scores_rows <- scores_rows[,-1]

#then sort within partition
i <- 2
partition <- 1
partitions <- list()
q <- 1
while (i <= ncol(currentmatrix)) {
  while (i <= ncol(currentmatrix) && isTRUE(all.equal(currentmatrix[,i],currentmatrix[,i-1]))) {
    partition <- c(partition,i)
    i <- i + 1
  }
  partitions[[q]] <-  partition
  partition <- i # set the new starting value for next partition
  i <- i+1 # set new i
  q <- q+1
}

if (is.null(dim(scores_rows)) | nrow(scores_rows) == 1) {
for (i in 1:length(partitions)) {
  if (length(partitions[[i]]) > 1) {
    scores_rows[,partitions[[i]]] <- t(apply(scores_rows[,partitions[[i]]],1,sort,decreasing = TRUE))
  }
}
}
else {
  for (i in 1:length(partitions)) {
    if (length(partitions[[i]]) > 1) {
      scores_rows[partitions[[i]]] <- sort(scores_rows[partitions[[i]]],decreasing = TRUE)
    }
  }
}
scores_rows <- vector.to.row(scores_rows) # making this a row again just in case
#then knock out the identical ones
newrows <- unique(split(as.matrix(scores_rows),row(scores_rows)))
newmatrices <- list()
for (i in 1:length(newrows)) {
  temp <- rbind(currentmatrix,newrows[[i]])
  temp <- temp[,(colSums(temp) != 0)]
  newmatrices[[i]] <- temp
}
newmatrices
}
