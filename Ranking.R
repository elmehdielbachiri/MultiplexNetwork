Ranking<-function(matrixrank,method){
  if(method=="bourda"){
    Bestindex=which.max(rowSums(matrixrank))
  } else if (method=="Parisa"){
    Bestindex=pari(matrixrank)
  } else if (method=="Rastin"){
    Bestindex=Rast(matrixrank)
  }
  return(Bestindex)}
