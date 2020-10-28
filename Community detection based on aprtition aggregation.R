multiplex.partition.aggregation <- function(multiplexObject,algorithm,h=10,alpha=0.7){
  Threshold <- alpha*h
  Threshold <- Threshold*multiplexObject@nbLayers
  mat <- matrix(c(0),ncol=length(multiplexObject@nodes),nrow=length(multiplexObject@nodes))
  res <- matrix(c(0),ncol=length(multiplexObject@nodes),nrow=length(multiplexObject@nodes))
  for(i in 1:h){
    for(graph in multiplexObject@layers){
      if(algorithm=="edge.betweenness"){
        wt <- edge.betweenness.community(graph)
        comList <- communities(wt)
      }
      
      else if(algorithm=="walktrap"){
        wt <- walktrap.community(graph)
        comList <- communities(wt)
      }
      
      else if(algorithm=="label.propagation"){
        wt <- label.propagation.community(graph)
        comList <- communities(wt)
      }
      
      else if(algorithm=="infomap"){
        wt <- infomap.community(graph)
        comList <- communities(wt)
      }
      
      else if(algorithm=="louvain"){
        wt <- multilevel.community(graph)
        comList <- communities(wt)
      }
      else if(algorithm=="licod"){
        comList <- licod(graph)
      }
      for(com in comList){
        for( x in 1:(length(com)-1)){
          for( y in (x+1):length(com)){
            mat[com[x],com[y]] <- (mat[com[x],com[y]] + 1)
            mat[com[y],com[x]] <- mat[com[x],com[y]]
          }
        }
      }
    }    
  }
  
  for(x in 1:(length(multiplexObject@nodes)-1)){
    for(y in (x+1):length(multiplexObject@nodes)){
      if( mat[x,y] >= Threshold){
        res[x,y] <- 1
        #res[x,y] <- mat[x,y]
      }
    }
  }
  
  graph_res <- graph.adjacency(res,mode = "undirected")
  clusters <- clusters(graph_res)
  return(communities(clusters))
}