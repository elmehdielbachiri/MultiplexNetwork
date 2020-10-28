rankList <- function(vect){
  vect <- rank(vect,ties.method = "min")
  tmp <- vector(mode="integer")
  for(i in 1:length(vect)){
    min <- min(vect[is.finite(vect)],na.rm=T)
    index <- which(vect==min)    
    if(length(index)==1){
      tmp[index] <- vect[index]
      vect[index] <- NA
    }
    else {
      for(j in 1:length(index)){
        tmp[index[j]] <- (vect[index[j]]+(j-1))
        vect[index[j]] <- NA
      }
    }
  }
  return(tmp)

}