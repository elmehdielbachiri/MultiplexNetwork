library(igraph)

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
heart_labels=heart[,2]

#Graph data set taht we convert to bidimensional matrices to create graphs from edgelists
zachary_dataset = read.csv('zachary',header = FALSE, skip=2,sep='',blank.lines.skip=TRUE)
zachary_matrix <- as.matrix(zachary_dataset)
zachary_graph <- graph_from_edgelist(zachary_matrix,directed=FALSE)

polbooks_dataset = read.csv('polbooks.mtx',header = FALSE, skip=2,sep='')
polbooks_matrix <- as.matrix(polbooks_dataset)
polbooks_graph <- graph_from_edgelist(polbooks_matrix,directed=FALSE)

polblogs_dataset = read.csv('polblogs.mtx',header = FALSE, skip=1,sep='')
polblogs_dataset$V3<-NULL
polblogs_matrix <- as.matrix(polblogs_dataset)
polblogs_graph <- graph_from_edgelist(polblogs_matrix,directed=FALSE)

dolphins_dataset = read.csv('soc-dolphins.mtx',header = FALSE, skip=1,sep='')
dolphins_matrix <- as.matrix(dolphins_dataset)
dolphins_graph <- graph_from_edgelist(dolphins_matrix,directed=FALSE)
