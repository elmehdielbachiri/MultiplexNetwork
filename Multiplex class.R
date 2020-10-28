library(igraph)
library(R.matlab)
library(base)

#'   Multiplex Class
setClass(
  # Set the name for the class
  "Multiplex",
  
  # Define the slots
  slots = c(
    nodes="vector", 
    nbLayers="numeric", 
    layers="list"
  ),
  
  # Set the default values for the slots
  prototype=list(
    nodes=0, 
    nbLayers=0, 
    layers=list()
  ),
  
  # Function to see if the multiplex is consistent
  
  validity = function(object){
    i = object@nbLayers
    while(i!=0){
      if(length(V(object@layers[[i]])$name)!=vcount(object@layers[[i]])){
        stop (" the number of vertex is not the same in all layers")
      }
    }
    return(TRUE)
  }
)
