library(igraph)
library(pipeR)

from <- sapply(1:31,function(x) rep(rownames(mtcars)[x],32-x)) %>>% unlist
to <- sapply(2:32,function(x) rownames(mtcars)[x:32]) %>>% unlist
d <- dist(mtcars[1:7])

q1 <- summary(d)[2]
# Make a graph/network
g <- graph_from_data_frame(data.frame(from=from[which(d <= q1)],
									  to=to[which(d <= q1)],
									  relations=as.vector(d)[which(d <= q1)]),
						   vertices=cbind(names=rownames(mtcars),mtcars[1:7]),
						   directed=F)

set.seed(1984)
louvain_cl <- cluster_louvain(g)
set.seed(1984)
fast_greedy_cl <- cluster_fast_greedy(g)
set.seed(1984)
edge_betweenness_cl <- cluster_edge_betweenness(g)
set.seed(1984)
infomap_cl <- cluster_infomap(g)
set.seed(1984)
label_prop_cl <- cluster_label_prop(g)
set.seed(1984)
leading_eigen_cl <- cluster_leading_eigen(g)
set.seed(1984)
walktrap_cl <- cluster_walktrap(g)
set.seed(1984)
pam_cl <- cluster::pam(mtcars[1:7],3)

library(cowplot)
library(ggsci)
set.seed(1984)
g1 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(louvain_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=15,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="louvain"))
set.seed(1984)
g2 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(fast_greedy_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=18,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="fast_greedy"))
set.seed(1984)
g3 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(edge_betweenness_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=19,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="edge_betweenness"))
set.seed(1984)
g4 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(infomap_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=15,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="infomap"))
set.seed(1984)
g5 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(label_prop_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=18,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="label_prop"))
set.seed(1984)
g6 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(leading_eigen_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=20,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="leading_eigen"))
set.seed(1984)
g7 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(walktrap_cl$membership));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=15,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="walktrap"))
set.seed(1984)
g8 <- (Rtsne::Rtsne(mtcars[1:7],perplexity=10)$Y) %>>% 
	{. <- data.frame(x=.[,1],
					 y=.[,2],
					 clusters=as.factor(pam_cl$cluster));.} %>>%
(ggplot(.,aes(x,y)) + geom_point(aes(color=clusters),shape=18,size=6,alpha=0.7) +
 scale_color_lancet() + theme_bw() + labs(title="PAM"))

mst <- minimum.spanning.tree(g)
library(ggraph)
g_mst <- ggraph(mst,'igraph',algorithm='tree',circular=T) + 
	geom_edge_diagonal(aes(alpha=..index..)) + coord_fixed() + 
	geom_node_point(aes(size=degree(mst)),color='steelblue') + 
	geom_node_text(aes(label=name),color='darkmagenta',size=3) + 
	ggforce::theme_no_axes()

show(plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g_mst))
