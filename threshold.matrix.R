threshold.matrix <- function(mat, threshold,weighted){
  A <- matrix(c(0), nrow=nrow(mat), ncol=ncol(mat))
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(mat[i,j]>=threshold){
        if(weighted)
              A[i,j] <- mat[i,j] 
        else 
          A[i,j] <- 1
      }
    }
  }
  return(A)
}

multiplex.threshold.matrix <- function(vect, threshold,weighted){
  if(! is.vector(vect,"list")){
    print("not a list ")
    exit(0)
  }  
  A <- vector("list",length(vect))
  for(i in 1:length(vect)){
  A[[i]] <- threshold.matrix(vect[[i]], threshold,weighted = weighted)
  }
  return(A)
}