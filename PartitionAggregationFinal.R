library(igraph)
library(R.matlab)
library(base)
library(igraph)
library(devtools)
load_all()
library(clue)
library(clusterSim)
library(clusterCrit)
library(cccd)
library(mclust)
library(sampling)
library(igraph)
library(clusteval)
library(graphics)

### Try normalizing and standardising data and quality indices, (udemy) or use scale
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
            Cons[is.na(Cons)] <- 0
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
multiplex.load.graph <- function(L,weighted = F){
  if(! is.list(L)){
    #print("not a list ")
    exit(0)
  }
  #create a multiplex object
  multiplexObject <- new("Multiplex")
  
  multiplexObject@nbLayers = length(L)  
  
  #read layers 
  for(i in 1:multiplexObject@nbLayers){
    ##print(get.adjacency(L[[i]]))
    multiplexObject@layers[[i]] <- graph.adjacency(adjmatrix = threshold.matrix(get.adjacency(L[[i]]),threshold = 1,weighted = weighted),mode = "undirected")
  }
  multiplexObject@nodes=c(1:length(V(multiplexObject@layers[[1]])))
  
  return(multiplexObject)
}



b<-iris[1:150,5]
a<-glass[1:150,11]


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

#####
Ranking<-function(matrixrank,method){
  if(method=="bourda"){
    Bestindex=which.max(rowSums(matrixrank))
  } else if (method=="Parisa"){
    Bestindex=pari(matrixrank)
  } else if (method=="Rastin"){
    Bestindex=Rast(matrixrank)
  }
  return(Bestindex)}

###### DATA

###### DATA


#Vector data set
iris = read.csv('iris.data',header=TRUE)
wine = read.csv('wine.data', header=TRUE)
glass= read.csv('glass.data', header = FALSE)
breastcancer = read.csv('wdbc.data', header = FALSE, sep=',')
soybean = read.csv('soybean-small.data',header=FALSE, sep=',')
heart = read.csv('processed.cleveland.data',header=FALSE, sep=',',na.strings = '?')
heart <- na.omit(heart)
# here we NULL out the column that gives the class (cluster membership)
iris_dataset = iris[1:4]
iris_labels = iris[,5]
wine_dataset = wine[2:14]
wine_labels = wine[,1]
glass_dataset = glass[1:10]
glass_labels = glass[,11]
breast_dataset = breastcancer[3:32]
breast_labels= breastcancer[,2]
soybean_dataset = soybean[1:35]
soybean_labels=soybean[,36]
heart_dataset=heart[-2]
heart_labels=heart[,14]

#Graph data set taht we convert to bidimensional matrices to create graphs from edgelists
karate <- as.undirected(read.graph("karate.gml", format="gml"))
karate_labels <- vertex_attr(karate,'value')

dolphins <- as.undirected(read.graph("dolphins.gml", format="gml"))
dolphins_labels <- vertex_attr(dolphins,'value')

polbooks <- as.undirected(read.graph("polbooks.gml", format="gml"))
polbooks_labels <- factor(vertex_attr(polbooks,'value'))

polblogs <- as.undirected(read.graph("polblogs.gml", format="gml"))
polblogs_labels <- vertex_attr(polblogs,'value')
#####Community detection
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

####### Consensus


## Librairies ##



FinalM2<-function(P=15,C=5,DATA,label,h=10,alpha=0.6,method="bourda"){
  # P = nombre de partitions
  
  # C = nombre de clusters ? trouver pour k-means
  # Data = les donn?es
  # label = les labels connus
  # method = methode ? utiliser pour le ranking
  
  
  ## GENERATION DES PARTITIONS INITIALES ##
  
  # P1 = une seule partition
  if(is_igraph(DATA)){
    P1 = t(membership(cluster_label_prop(DATA))) # si c'est un graphe, applique Label Propagation
  }
  else{
    P1 = kmeans(DATA,C)$cluster # si c'est des données vectorielles, applique Kmeans
  }
  
  # M = Matrice qui contient toutes les partitions de l'ensemble initial
  M=matrix(0,length(label),P)
  for (n in 1:P){
    if(is_igraph(DATA)){
      M[,n]<-t(membership(cluster_label_prop(DATA))) # si c'est un graphe, applique Label Propagation
    }
    else{
      # Matrice de PartitioN
      M[,n]<-kmeans(DATA,C)$cluster  ## si c'est des données vectorielles, applique Kmeans
    }
  }
  
  ## CREATION DU MULTIPLEX AVEC DIFFERENTS INDICES DE DIVERSITE ##
  
  ####Normalize "split.join" 
  L<-Multiplex(M,c("VI","rand","adjusted.rand","nmi","jaccard")) #list des graph du multiplex
  Multi<-multiplex.load.graph(L)# creation du multiplex
  #Multi = graph multiplexe
  
  ## DETECTION DES COMMUNAUTES DANS LE MULTIPLEX ##
  
  # Comm = les communautés dans le mutiplex
  #Comm<-muxLicod(Multi) # the communities detection (la methode Muxlicod)
  # Genlouvain 
  Comm<-multiplex.partition.aggregation(Multi,algorithm ="louvain",h,alpha) # the communities detection (la method Louvain)
  
  nPartitionApresSelection=(length(Comm))
  
  
  
  ## CALCUL DES INDICES DE QUALITE ##
  
  # Qtot = la matrice des indices de qualité
  Qtot=matrix(0,P,5)
  
  for(k in 1:P){
    l=1
    for(s in Comm){ # trouve les label de partitions, elles sont dans quelle communauté
      if(k %in% s){
        r=l
      }
      l=l+1
    }
    
    Qtot[k,1]=k # numero de partition
    Qtot[k,2]=r # label de partition, elle est dans quelle communite
    if(is_igraph(DATA)){
      Qtot[k,3]= modularity(DATA, M[,k]) #applique indexe de quality Modularity
      Qtot[k,4]= modularity(DATA, M[,k]) #applique indexe de quality Modularity
      Qtot[k,5]= modularity(DATA, M[,k]) #applique indexe de quality Modularity
    }
    else{
      Qtot[k,3]=1/(index.DB(DATA,M[,k])$DB) #applique indexe de quality Davies-Bouldin
      Qtot[k,4]=index.G1(DATA,M[,k]) #applique indexe de quality Calinski-Harabasz
      d <- dist(DATA)
      Qtot[k,5]=(index.S(d,M[,k])+1)/2 #applique indexe de quality Sillouette
    }
  }
  for (i in 3:5){
    Qtot[,i] <- Qtot[,i]/max(Qtot[,i])
  }
  #print(Qtot)
  
  
  ## SELECTION FINALE AVEC LES INDICES DE QUALITE (RANK) ##
  
  
  # FinalPartitionIndex = les numero des partitions de meilleure qualité selon le ranking
  FinalPartitionIndex=matrix(0,length(Comm),1)
  for(d in 1:length(Comm)){ #construit une matrix de Ranking
    G=Qtot[Qtot[,2]==d,]
    if(is.null(dim(G)[1])){
      FinalPartitionIndex[d,1]=G[1]
    }
    else{
      matrixrank=G[,1:3]
      matrixrank[,1]=rank(G[,3])
      matrixrank[,2]=rank(G[,4])
      matrixrank[,3]=rank(G[,5])
      FinalPartitionIndex[d,1]=G[Ranking(matrixrank,method),1] # Fait le ranking
    }
  }
  # FinalPartition = matrice des partitions apres selection
  FinalPartition=M[,FinalPartitionIndex]
  
  
  ## CONSENSUS APRES SELECTION ##
  
  # PartitionConsensus = la partion finale apres consensus
  
  dimp=dim(FinalPartition)[2] # dimp = nombre de partitions selectionées
  
  if (is.null(dimp)){ # si une seule partiton selectionée, pas de consensus
    PartitionConsensus = FinalPartition
  }
  
  else{ # sinon on fait le consensus
    matrixData=matrix(0,length(label),length(label)) # matrice du graphe (length(label) = nombre de donn?es)  
    for(h in 1:dimp){ # pour chaque partition selectionée
      for(i in 1:(length(label)-1)){ # pour chaque paire de données
        for(j in (i+1):length(label)){
          if (FinalPartition[i,h]==FinalPartition[j,h]){ # si les deux donn?es sont dans le meme cluster
            matrixData[i,j]=matrixData[i,j]+1 # on renforce le lien entre ces deux donn?es dans le graphe
            matrixData[j,i]=matrixData[i,j]
          }
        }
      }
    }
    
    matrixData=graph.adjacency(matrixData, mode="undirected",weighted=TRUE) # construction du graphe de consensus
    PartitionConsensus=walktrap.community(matrixData)$membership # detection de community dans le graphe pour obtenir la partition finale
  }
  
  
  
  ## CONSENSUS SANS SELECTION ##
  
  # pareil qu'avant mais avec toutes les partitions
  matrixDataSansSelection=matrix(0,length(label),length(label))
  for(h in 1:dim(M)[2]){
    for(i in 1:(length(label)-1)){ # pour chaque paire de données
      for(j in (i+1):length(label)){
        if (M[i,h]==M[j,h]){
          matrixDataSansSelection[i,j]=matrixDataSansSelection[i,j]+1
          matrixDataSansSelection[j,i]=matrixDataSansSelection[i,j]
        }
      }
    }
  }
  
  matrixDataSansSelection=graph.adjacency(matrixDataSansSelection, mode="undirected",weighted=TRUE) # construction du graphe de consensus
  PartitionConsensusSansSelection=walktrap.community(matrixDataSansSelection)$membership # detection de community dans le graphe pour obtenir la partition finale
  
  
  
  ## VALEURS A RETOURNER (Adj RAND, NMI, ETC) ##
  
  NMInos = compare(t(PartitionConsensusSansSelection),label,method="nmi")
  NMI = compare(t(PartitionConsensus),label,method="nmi")
  ARInos = compare(t(PartitionConsensusSansSelection),label,method="adjusted.rand")
  ARI = compare(t(PartitionConsensus),label,method="adjusted.rand")
  NMIP1 = compare(t(P1),label,method="nmi")
  ARIP1 = compare(t(P1),label,method="adjusted.rand")
  nclustP1=max(P1)
  nclustPC=max(PartitionConsensus)
  nclustPCSS=max(PartitionConsensusSansSelection)
  
  return(list(P1 = P1, PC = PartitionConsensus,PCSS =PartitionConsensusSansSelection,ARSS=ARInos,ARI=ARI, ARIP1 = ARIP1,NMISS=NMInos,NMI=NMI,NMIP1 =NMIP1,nclustP1=nclustP1,dimp,nclustPC=nclustPC,nclustPCSS=nclustPCSS,nPartitionApresSelection=nPartitionApresSelection))
}

Y=FinalM2(P=10,C=5,DATA=soybean_dataset,soybean_labels,10,0.7); Y



finaltest <- function(nbofpartitions,nbofexecution=35,D,lab,centers=5){
  SPNMI <- c()
  ECNMI <- c()
  CESNMI <- c()
  SPARI <- c()
  ECARI <- c()
  CESARI <- c()
  for(i in 1:nbofexecution){
    #for(k in nbpart){
    Y=FinalM2(P=nbofpartitions,C=centers,D,lab,10,0.7)
    #Simple partitioning NMI
    SPNMI <- c(SPNMI,Y$NMIP1)
    ECNMI <- c(ECNMI,Y$NMISS)
    CESNMI <- c(CESNMI,Y$NMI)
    SPARI <- c(SPARI,Y$ARIP1)
    ECARI <- c(ECARI,Y$ARSS)
    CESARI <- c(CESARI,Y$ARI)
    #}
  } 
  moyenneSPNMI <- mean(SPNMI,na.rm = TRUE)
  moyenneECNMI <- mean(ECNMI,na.rm = TRUE)
  moyenneCESNMI <- mean(CESNMI,na.rm = TRUE)
  moyenneSPARI <- mean(SPARI,na.rm = TRUE)
  moyenneECARI <- mean(ECARI,na.rm = TRUE)
  moyenneCESARI <- mean(CESARI,na.rm = TRUE)
  
  ectypeSPNMI <- sd(SPNMI,na.rm = TRUE)
  ectypeECNMI <- sd(ECNMI,na.rm = TRUE)
  ectypeCESNMI <- sd(CESNMI,na.rm = TRUE)
  ectypeSPARI <- sd(SPARI,na.rm = TRUE)
  ectypeECARI <- sd(ECARI,na.rm = TRUE)
  ectypeCESARI <- sd(CESARI,na.rm = TRUE)
  return(list(c(moyenneSPNMI,ectypeSPNMI),c(moyenneECNMI,ectypeECNMI),c(moyenneCESNMI,ectypeCESNMI),c(moyenneSPARI,ectypeSPARI),c(moyenneECARI,ectypeECARI),c(moyenneCESARI,ectypeCESARI)))
}


#finaltest(35,nbofexecution=30,dolphins,dolphins_labels,centers=2)


alldatatest <- function (nbofexecution=35,nbofpartitions,centers){
  # 6 datasets
  A=matrix(0,length(VECTDAT),13)
  for (i in 1:length(VECTDAT)){
    # Try this as well P=15, C=4, sigma=0.5, delta=0.5, theta=0.6, iter=1
    Y=finaltest(nbofpartitions,nbofexecution,VECTDAT[[i]],VECTLABELS[[i]],centers[[i]])
    A[i,1] <- i
    A[i,2] <- Y[[1]][1]
    A[i,3] <- Y[[1]][2]
    A[i,4] <- Y[[2]][1]
    A[i,5] <- Y[[2]][2]
    A[i,6] <- Y[[3]][1]
    A[i,7] <- Y[[3]][2] 
    A[i,8] <- Y[[4]][1]
    A[i,9] <- Y[[4]][2]
    A[i,10] <- Y[[5]][1]
    A[i,11] <- Y[[5]][2]
    A[i,12] <- Y[[6]][1]
    A[i,13] <- Y[[6]][2]
  }
  return(A)
}


#alldatatest(nbofexecution=20,15,centers=centersvect)




finaltest <- function(nbofpartitions,nbofexecution=35,D,lab,centers=5){
  SPNMI <- c()
  ECNMI <- c()
  CESNMI <- c()
  SPARI <- c()
  ECARI <- c()
  CESARI <- c()
  for(i in 1:nbofexecution){
    #for(k in nbpart){
    Y=FinalM2(P=nbofpartitions,C=centers,D,lab,30,0.7)
    #Simple partitioning NMI
    SPNMI <- c(SPNMI,Y$NMIP1)
    ECNMI <- c(ECNMI,Y$NMISS)
    CESNMI <- c(CESNMI,Y$NMI)
    SPARI <- c(SPARI,Y$ARIP1)
    ECARI <- c(ECARI,Y$ARSS)
    CESARI <- c(CESARI,Y$ARI)
    #}
  } 
  moyenneSPNMI <- mean(SPNMI)
  moyenneECNMI <- mean(ECNMI)
  moyenneCESNMI <- mean(CESNMI)
  moyenneSPARI <- mean(SPARI)
  moyenneECARI <- mean(ECARI)
  moyenneCESARI <- mean(CESARI)
  
  ectypeSPNMI <- sd(SPNMI,na.rm = TRUE)
  ectypeECNMI <- sd(ECNMI,na.rm = TRUE)
  ectypeCESNMI <- sd(CESNMI,na.rm = TRUE)
  ectypeSPARI <- sd(SPARI,na.rm = TRUE)
  ectypeECARI <- sd(ECARI,na.rm = TRUE)
  ectypeCESARI <- sd(CESARI,na.rm = TRUE)
  return(list(c(moyenneSPNMI,ectypeSPNMI),c(moyenneECNMI,ectypeECNMI),c(moyenneCESNMI,ectypeCESNMI),c(moyenneSPARI,ectypeSPARI),c(moyenneECARI,ectypeECARI),c(moyenneCESARI,ectypeCESARI)))
}



VECTDAT=list(heart_dataset)
VECTLABELS=list(heart_labels)
centersvect=list(4)
  
  
alldatatest <- function (nbofexecution=35,nbofpartitions,centersvect){
  # 6 datasets
  A=matrix(0,length(VECTDAT),13)
  for (i in 1:length(VECTDAT)){
    # Try this as well P=15, C=4, sigma=0.5, delta=0.5, theta=0.6, iter=1
    Y=finaltest(nbofpartitions,nbofexecution=5,VECTDAT[[i]],VECTLABELS[[i]],centersvect[[i]])
    A[i,1] <- i
    A[i,2] <- Y[[1]][1]
    A[i,3] <- Y[[1]][2]
    A[i,4] <- Y[[2]][1]
    A[i,5] <- Y[[2]][2]
    A[i,6] <- Y[[3]][1]
    A[i,7] <- Y[[3]][2] 
    A[i,8] <- Y[[4]][1]
    A[i,9] <- Y[[4]][2]
    A[i,10] <- Y[[5]][1]
    A[i,11] <- Y[[5]][2]
    A[i,12] <- Y[[6]][1]
    A[i,13] <- Y[[6]][2]
  }
  return(A)
}


alldatatest(nbofexecution=35,30,centers=centersvect)




