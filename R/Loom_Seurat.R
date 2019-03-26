#library(Seurat)
library(ggplot2)
source("./R/Seurat_functions.R")
library(Seurat)
#====== 6.1 load seurat  ==========================================
lnames = load(file = "./data/DLBCL_8_alignment.Rda")
head(DLBCL@meta.data)
table(DLBCL@meta.data$orig.ident)
lnames
DLBCL.subsets <- SplitCells(DLBCL)
DLBCL.subsets[[3]]
#DLBCL.PDX <- DLBCL.subsets[[1]]
DLBCL.primary <- DLBCL.subsets[[2]]
remove(DLBCL)
remove(DLBCL.subsets)
GC()
idents <- as.data.frame(table(DLBCL.primary@ident))
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
DLBCL.primary@ident <- plyr::mapvalues(x = DLBCL.primary@ident,
                                   from = old.ident.ids,
                                   to = new.cluster.ids)

gg_colors <- function(object = DLBCL.primary, no.legend = TRUE, do.label = TRUE,
                      do.return = TRUE, label.size = 6, gg_title="primary DLBCL"){
        g1 <- Seurat::TSNEPlot(object = object, no.legend = no.legend,
                               do.label = do.label,do.return = do.return,
                               label.size = label.size)+
                ggtitle(gg_title)+
                theme(text = element_text(size=15),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
        print(g1)
        g <- ggplot2::ggplot_build(g1)
        print(unique(g$data[[1]]["colour"]))
        colors <- unlist(g$data[[1]]["colour"])
        cells <- Seurat::WhichCells(object)
        names(colors) <- stringi::stri_replace_last_fixed(cells,"_",":")
        names(colors) <- paste0(names(colors),"x")
        
        return(colors)
}
primary_colors <- gg_colors(object = DLBCL.primary)
#PDX_colors <- gg_colors(object = DLBCL.PDX)

#select TSNE embeddings
DLBCL.primary_emb <- DLBCL.primary@dr$tsne@cell.embeddings[,1:2]
row.names(DLBCL.primary_emb) <- stringi::stri_replace_last_fixed(row.names(DLBCL.primary_emb),
                                                          "_",":")
row.names(DLBCL.primary_emb) <- paste0(row.names(DLBCL.primary_emb),"x")
remove(DLBCL.primary)

#DLBCL.PDX_emb <- DLBCL.PDX@dr$tsne@cell.embeddings[,1:2]
#row.names(DLBCL.PDX_emb) <- stringi::stri_replace_last_fixed(row.names(DLBCL.PDX_emb),
#                                                                 "_",":")
#row.names(DLBCL.PDX_emb) <- paste0(row.names(DLBCL.PDX_emb),"x")
#remove(DLBCL.PDX)

#====== 6.2 load primary loom  ==========================================
ldat  <- velocyto.R::read.loom.matrices("./data/primary_3.loom")
colors <- primary_colors
emb <- DLBCL.primary_emb

#PDX_dat  <- velocyto.R::read.loom.matrices("./data/PDX_5.loom")
#colors <- PDX_colors
#emb <- DLBCL.PDX_emb
save(ldat,colors,emb, file = "./data/primary_3.Rda")
GC()
################################################################
# Gene filtering
# Spliced expression magnitude distribution across genes:
#####################################################
lnames = load(file = "./data/primary_3.Rda")
lnames
library(velocyto.R)
hist(log10(Matrix::rowSums(ldat$spliced)+1),col='wheat',
     xlab='log10[ number of reads + 1]',
     main='number of reads per gene')
#Set up expression matrices, filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$ambiguous
# filter expression matrices based on some minimum max-cluster averages

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,colors,min.max.cluster.average = 0.1)
nmat <- filter.genes.by.cluster.expression(nmat,colors,min.max.cluster.average = 0.1)
smat <- filter.genes.by.cluster.expression(smat,colors,min.max.cluster.average = 0.1)
# look at the resulting gene set
str(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

# We calculate gene-relative velocity, 
# using k=5 cell kNN pooling, but now using entire range of expression to determine slope gamma,
# and using spanning reads (smat) to fit the gene offsets.
rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02)
pca.velocity.plot(rvel.qf,nPcs=5,plot.cols=2,cell.colors=ac(colors,alpha=0.7),
                  cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1))


rvel1 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,deltaT2 = 1,kCells = 1)
pca.velocity.plot(rvel1,nPcs=5,plot.cols=2,cell.colors=ac(colors,alpha=0.7),
                  cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,1,1,1))

#Visualization on an existing embedding
rvel <- rvel.qf;arrow.scale=3; cell.alpha=0.4; cell.cex=1
show.velocity.on.embedding.cor(emb,rvel,n=100,scale='sqrt',
                               cell.colors=ac(colors,alpha=cell.alpha),
                               cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)
arrow.scale=8
show.velocity.on.embedding.cor(emb,rvel,n=100,scale='sqrt',
                               cell.colors=ac(colors,alpha=cell.alpha),
                               cex=cell.cex,arrow.scale=arrow.scale,
                               show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=20,arrow.lwd=2)
