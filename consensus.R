## Librairies ##

library(igraph)
library(devtools)
load_all()
library(clue)
library(clusterSim)
library(clusterCrit)
library(cccd)
library(mclust)
#source('C:/Users/paris/Dropbox/Travail_2018_19/Article/CES_Jornal_Paper/Matrix diss.R') #matrix de dissimilarité
#source('C:/Users/paris/Dropbox/Travail_2018_19/Article/CES_Jornal_Paper/rankinglist.R') #Matrix de ranking

## Fonction pour les experiences ##

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
  print(Qtot)
  
  
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

Y=FinalM2(P=100,C=5,DATA=soybean_dataset,soybean_labels,10,0.6); Y


#Criteria using clustersim
intCriteria(cbind(seq(length(iris_labels)),Y$PC),c(iris_labels),'all')

VECTDAT=list(iris_dataset,
             wine_dataset,
             glass_dataset,
             breast_dataset,
             soybean_dataset,
             heart_dataset)
VECTLABELS=list(iris_labels,
                wine_labels,
                glass_labels,
                breast_labels,
                soybean_labels,
                heart_labels)

testpartitionaggregation <- function(n){
  ariCOUNTERlist <- c()
  nmiCOUNTERlist <- c()
  for (i in 1:length(VECTDAT)){
    ariCOUNTER=0
    nmiCOUNTER=0
    for (k in 1:n){
      Y=FinalM2(P=15,C=5,DATA=VECTDAT[[i]],VECTLABELS[[i]],20,0.7)
      if(Y$ARI >= Y$ARSS & Y$ARI >= Y$ARIP1){
        ariCOUNTER <- ariCOUNTER+1}
      if(Y$NMI >= Y$NMISS & Y$NMI >= Y$NMIP1){
        nmiCOUNTER <- nmiCOUNTER+1}
    }
    ariCOUNTERlist<-c(ariCOUNTERlist,ariCOUNTER) 
    nmiCOUNTERlist<-c(nmiCOUNTERlist,nmiCOUNTER) 
  }  
  return(list(ariCOUNTERlist,nmiCOUNTERlist))
}





# 6 datasets
A=matrix(0,length(VECTDAT),7)
for (i in 1:length(VECTDAT)){
  Y=FinalM2(P=15,C=5,DATA=VECTDAT[[i]],VECTLABELS[[i]],100,0.7)
  A[i,1] <- i
  A[i,2] <- Y$NMIP1
  A[i,3] <- Y$ARIP1
  A[i,4] <- Y$NMISS
  A[i,5] <- Y$ARSS
  A[i,6] <- Y$NMI
  A[i,7] <- Y$ARI
  print(A[i,])
}


# influence of alpha
indexP1=c()
indexPC=c()
indexPCSS=c()
for (h in seq(1,50)){
  Y <- FinalM2(P=15,C=5,DATA=soybean_dataset,soybean_labels,h,alpha=0.7)
  indexP1 <- c(indexP1,Y$NMIP1)
  indexPCSS <- c(indexPCSS,Y$NMISS)
  indexPC <- c(indexPC,Y$NMI)
}
plot(seq(1,50),indexP1,col="black",pch="o",lty=1,ylim=c(0.4,1.2),main="influence of iterations h",ylab ="NMI index soybean")
lines(seq(1,50), indexP1, col="black",lty=1)
points(seq(1,50), indexPCSS, col="blue", pch="*")
lines(seq(1,50), indexPCSS, col="blue",lty=2)
points(seq(1,50), indexPC, col="dark red",pch="+")
lines(seq(1,50), indexPC, col="dark red", lty=3)
legend(0,1.2,legend=c("Partition Consensus","Partition Consensus with selection","Simple partition"), col=c("blue","red","black"),
       pch=c("o","*","+"),lty=c(1,2,3), ncol=1)


# Chercher labels
