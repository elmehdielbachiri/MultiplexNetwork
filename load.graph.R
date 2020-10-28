multiplex.load.graph <- function(L,weighted = F){
            if(! is.list(L)){
              print("not a list ")
              exit(0)
            }
            #create a multiplex object
            multiplexObject <- new("Multiplex")
            
            multiplexObject@nbLayers = length(L)  
            
            #read layers 
            for(i in 1:multiplexObject@nbLayers){
              #print(get.adjacency(L[[i]]))
              multiplexObject@layers[[i]] <- graph.adjacency(adjmatrix = threshold.matrix(get.adjacency(L[[i]]),threshold = 1,weighted = weighted),mode = "undirected")
            }
            multiplexObject@nodes=c(1:length(V(multiplexObject@layers[[1]])))
            
            return(multiplexObject)
          }


