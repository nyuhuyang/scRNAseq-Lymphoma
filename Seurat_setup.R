########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
# remove "hg19_" tag from genes.tsv files
# open shell
# cd ./data/IP8672/outs/filtered_gene_bc_matrices/hg19
# sed 's/hg19_//g' genes.tsv > genes1.tsv
# mv genes1.tsv genes.tsv
# cd ./data/TR624_M2102/outs/filtered_gene_bc_matrices/hg19
# sed 's/hg19_//g' genes.tsv > genes1.tsv
# mv genes1.tsv genes.tsv

#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

samples <- c("IP8672","TR624_M2102")
projects <- c("EC-RW-4262","EC-RW-4262")
conditions <- c("primary", "PDX")
DLBCL_raw <- list()
for(i in 1:length(samples)){
  DLBCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                              samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
  colnames(DLBCL_raw[[i]]) <- paste0(conditions[i],
                                          "_",colnames(DLBCL_raw[[i]]))
}
DLBCL_Seurat <- lapply(DLBCL_raw, CreateSeuratObject,
                            min.cells = 3,
                            min.genes = 200,
                            project = projects)
for(i in 1:length(samples)) DLBCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
DLBCL_Seurat <- lapply(DLBCL_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 500, 
                            high.thresholds = Inf)
DLBCL_Seurat <- lapply(DLBCL_Seurat, NormalizeData)
DLBCL_Seurat <- lapply(DLBCL_Seurat, ScaleData)
DLBCL_Seurat <- lapply(DLBCL_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(DLBCL_Seurat, function(x) head(rownames(x@hvg.info), 1000))
genes.use <- unique(c(g[[1]],g[[2]]))
for(i in 1:length(conditions)){
  genes.use <- intersect(genes.use, rownames(DLBCL_Seurat[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size 16764

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
DLBCL <- RunCCA(DLBCL_Seurat[[1]],DLBCL_Seurat[[2]],
                          genes.use = genes.use,
                          num.cc = 30)
save(DLBCL, file = "./data/DLBCL_alignment.Rda")
remove(DLBCL_raw)
remove(DLBCL_Seurat)

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = DLBCL, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = DLBCL, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = DLBCL, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

p3 <- MetageneBicorPlot(DLBCL, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE)
p3 + geom_smooth(method = 'loess')
DimHeatmap(object = DLBCL, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)

DimHeatmap(object = DLBCL, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

PrintDim(object = DLBCL, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)


#======1.3 QC (skip)==================================


#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

DLBCL <- AlignSubspace(object = DLBCL, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:20)
#Now we can run a single integrated analysis on all cells!
DLBCL <- RunTSNE(object = DLBCL, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
DLBCL <- FindClusters(object = DLBCL, reduction.type = "cca.aligned", dims.use = 1:20, 
                      resolution = 0.8, force.recalc = T, save.SNN = TRUE)
p1 <- TSNEPlot(DLBCL, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(DLBCL, do.label = F, do.return = T, pt.size = 1)

plot_grid(p1, p2)


#Now, we annotate the clusters as before based on canonical markers.

TSNEPlot(object = DLBCL,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
  ggtitle("TSNE plot of all clusters")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle

# Compare clusters for each dataset
cell.all <- FetchData(DLBCL,"conditions")
cell.primary <- rownames(cell.all)[cell.all$conditions =="primary"]
cell.PDX <- rownames(cell.all)[cell.all$conditions =="PDX"]

DLBCL.primary <- SubsetData(object = DLBCL,
                               cells.use =cell.primary)
DLBCL.PDX <- SubsetData(object = DLBCL,
                              cells.use =cell.PDX)
table(DLBCL.primary@ident)
table(DLBCL.PDX@ident)
p1 <- TSNEPlot(object = DLBCL.primary,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
    ggtitle("Primary sample")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle

p2 <- TSNEPlot(object = DLBCL.PDX,do.label = TRUE, group.by = "ident", 
               do.return = TRUE, no.legend = TRUE,
               pt.size = 1,label.size = 8 )+
    ggtitle("PDX sample")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)
save(DLBCL, file = "./data/DLBCL_alignment.Rda")
remove(DLBCL.PDX)
remove(DLBCL.primary)
