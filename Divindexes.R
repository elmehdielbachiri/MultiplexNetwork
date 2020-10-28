library(igraph)
library(clusteval)

b<-iris[1:150,5]
a<-glass[1:150,11]

# Normalized Mutual Information (NMI) measure
compare(a, b, method = c("nmi"))

# Variation of Information (VI) metric 2003:
compare(a, b, method = c("vi")) 

# Jaccard Index 2002:
clusteval::cluster_similarity(a, b, similarity = c("jaccard"), method = "independence") 

# van Dongen S metric 2000:
compare(a, b, method = c("split.join")) 

# Adjusted Rand Index 1985:
compare(a, b, method = c("adjusted.rand")) 

# Rand Index 1971:
compare(a, b, method = c("rand")) 

# Jaccard Index:
compare2 <- function(c1,c2,method){
  if(method=="jaccard"){
  cluster_similarity(c1,c2,similarity = c("jaccard"))
  }
  else if(method=="vi"){
  1-compare(c1,c2,method = c("vi")) 
  }
  else{
    compare(c1,c2,method=c(method))
  }
  }
# Define normalized split.join :
normalized.split.join <- function(a,b){
  (split.join(a,b)+split.join(b,a))/(2*max(split.join(a,b)+split.join(b,a)))
}