library(sampling)
library(igraph)
library(clusteval)
library(graphics)

## creation des matrix de dissimilarité pour creer le multiplex avec indice de Diver ##

Multiplex<-function(M,Ind){
  
  # M = Matrice qui contient toutes les partitions de l'ensemble initial
  # Ind = c("Rand", "Adjusted.rand" ,"VI","split.join","NMI") : les indices ? utiliser
  
  P = dim(M)[2] # nombre de partitions
  K=length(Ind) #nmbr de indice de diversit?
  Cons=array(1,dim=c(P,P,K)) # matrice des diversité entre chaque partitions pour chaque indices de diversité
  for (i in 1:(P-1)){
    for (j in (i+1):P){
      for(k in 1:K){ # calcul des diversité entre chaque paire de partition et pour chaque indices
        Cons[i,j,k]=compare2(M[,i], M[,j], method = c(Ind[k])) 
        Cons[j,i,k]=Cons[i,j,k]
      }
    }
    
    # construction du Graph proximity basé sur relative neighbourhood graph
    
    G=array(1,dim=c(P,P,K))
    for (i in 1:(P-1)){
      for (j in (i+1):P){
        for (m in 1:P){
          for (k in 1:K){
            d<-Cons[i,j,k]
            dmi<-Cons[i,m,k]
            dmj<-Cons[m,j,k]
            if(d<dmi && d<dmj && m!=i && m!=j){
              G[i,j,k]=0
              G[j,i,k]=0
            }
          }
        }
      }
    }
  }
  
  for(i in 1:P){
    G[i,i,]=0 # pour eviter la boucle sur chaque node
  }
  
  # Creation du multiplex à partir de l'ensemble des graphes et affichage du multiplex
  
  g=list()
  for(k in 1:K){
    g[[k]] <- graph.adjacency(G[,,k], mode="undirected")
    #plot(g[[k]]) # affichage
    #title(main = Ind[k])
  }
  return(g)
}
