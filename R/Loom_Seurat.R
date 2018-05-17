library(velocyto.R)
library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 6.1 load seurat  ==========================================
lnames = load(file = "./data/DLBCL_8.Rda")
lnames
AllCells <- WhichCells(DLBCL)
IP8672.cells <- AllCells[grepl("IP8672",AllCells)]
IP8672 <- SubsetData(DLBCL, cells.use = IP8672.cells)
idents <- as.data.frame(table(IP8672@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "B cells",
                     "T & NK cells",
                     "B cells",
                     "Myeloid\n cells",
                     "T & NK cells",
                     "B cells")
IP8672@ident <- plyr::mapvalues(x = IP8672@ident,
                                   from = old.ident.ids,
                                   to = new.cluster.ids)

g1 <- TSNEPlot(object = IP8672, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("IP8672")+
        theme(text = element_text(size=15),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
print(g1)
g <- ggplot_build(g1)
unique(g$data[[1]]["colour"])
IP8672_colors <- unlist(g$data[[1]]["colour"])
names(IP8672_colors) <- gsub("_",":",IP8672.cells)
names(IP8672_colors) <- paste0(names(IP8672_colors),"x")
remove(DLBCL)
#select TSNE embeddings
IP8672_emb <- IP8672@dr$tsne@cell.embeddings[,1:2]
row.names(IP8672_emb) <- gsub("_",":",row.names(IP8672_emb))
row.names(IP8672_emb) <- paste0(row.names(IP8672_emb),"x")
#====== 6.2 load loom  ==========================================
ldat  <- read.loom.matrices("./data/IP8672/velocyto/IP8672.loom")

################################################################
# Gene filtering
# Spliced expression magnitude distribution across genes:
hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
#Set up expression matrices, filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$ambiguous
# filter expression matrices based on some minimum max-cluster averages

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,IP8672_colors,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,IP8672_colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,IP8672_colors,min.max.cluster.average = 0.5)
# look at the resulting gene set
str(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

# We calculate gene-relative velocity, 
# using k=5 cell kNN pooling, but now using entire range of expression to determine slope gamma,
# and using spanning reads (smat) to fit the gene offsets.
rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02)
pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cell.colors=ac(IP8672_colors,alpha=0.7),
                  cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1))


rvel1 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,deltaT2 = 1,kCells = 1)
pca.velocity.plot(rvel1,nPcs=5,plot.cols=2,cell.colors=ac(IP8672_colors,alpha=0.7),
                  cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,1,1,1))

#Visualization on an existing embedding
rvel <- rvel.qf;arrow.scale=3; cell.alpha=0.4; cell.cex=1
show.velocity.on.embedding.cor(IP8672_emb,rvel,n=100,scale='sqrt',
                               cell.colors=ac(IP8672_colors,alpha=cell.alpha),
                               cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)
arrow.scale=8
show.velocity.on.embedding.cor(IP8672_emb,rvel,n=100,scale='sqrt',
                               cell.colors=ac(IP8672_colors,alpha=cell.alpha),
                               cex=cell.cex,arrow.scale=arrow.scale,
                               show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=20,arrow.lwd=2)
