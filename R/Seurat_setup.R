########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(cowplot)

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
# Load the DLBCL dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

<<<<<<< HEAD

=======
<<<<<<< HEAD:R/Seurat_setup.R
>>>>>>> be7e28f4b38f6d7e5659c2f3478410bf014c01f3
samples <- c("3119T1", "3119T6","I17_1733", "IP-10927",
             "IP8672", "TR624_M2102", "TR_619_mouse_2080", "TR_619_mouse_2083")
projects <- paste0("EC-RW-",c(rep(4311,4),4262,4262,4311,4311))
conditions <- c("PDX", "PDX","primary","primary",
<<<<<<< HEAD
=======
=======
samples <- c("3119PR", "3119T1", "3119T6", "3139PR", "I17_1733", "IP-10927",
             "IP8672", "TR624_M2102", "TR_619_mouse_2080", "TR_619_mouse_2083")
projects <- paste0("EC-RW-",c(rep(4311,6),4262,4262,4311,4311))
conditions <- c("primary", "PDX", "PDX","primary","primary","primary",
>>>>>>> 9b8169132a35923308ef95843457f2dbd21713a4:R/Seurat_setup.R
>>>>>>> be7e28f4b38f6d7e5659c2f3478410bf014c01f3
                "primary","PDX","PDX","PDX")
DLBCL_raw <- list()
for(i in 1:length(samples)){
  DLBCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                              samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
  colnames(DLBCL_raw[[i]]) <- paste0(samples[i],"_",colnames(DLBCL_raw[[i]])) # Must be samples!
  rownames(DLBCL_raw[[i]]) <- sub("hg19_","",rownames(DLBCL_raw[[i]]))
}
DLBCL_Seurat <- lapply(DLBCL_raw, CreateSeuratObject,
                            min.cells = 2,
                            min.genes = 100,
                            project = projects)
for(i in 1:length(samples)) DLBCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
DLBCL_Seurat <- lapply(DLBCL_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 100, 
                            high.thresholds = Inf)
DLBCL_Seurat <- lapply(DLBCL_Seurat, NormalizeData)
DLBCL_Seurat <- lapply(DLBCL_Seurat, ScaleData)
DLBCL_Seurat <- lapply(DLBCL_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(DLBCL_Seurat, function(x) head(rownames(x@hvg.info), 500))
genes.use <- unique(unlist(g))
for(i in 1:length(conditions)){
  genes.use <- intersect(genes.use, rownames(DLBCL_Seurat[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size 16764

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
remove(DLBCL_raw)
WGCNA::collectGarbage()
DLBCL <- RunMultiCCA(object.list = DLBCL_Seurat, 
                          genes.use = genes.use,
                     niter = 25, num.ccs = 30,
                     standardize =TRUE)
<<<<<<< HEAD

=======
<<<<<<< HEAD:R/Seurat_setup.R
>>>>>>> be7e28f4b38f6d7e5659c2f3478410bf014c01f3
#save(DLBCL,genes.use, file = "./data/DLBCL_10_alignment.Rda")
DLBCL_marged_mat <- as(object = DLBCL@raw.data, "sparseMatrix")
DLBCL_metadata <- DLBCL@meta.data
save(DLBCL_marged_mat,DLBCL_metadata,genes.use, file = "./data/DLBCL_8.Rda")
<<<<<<< HEAD
=======

=======
save(DLBCL, file = "./data/DLBCL_10_alignment.Rda")
>>>>>>> 9b8169132a35923308ef95843457f2dbd21713a4:R/Seurat_setup.R
>>>>>>> be7e28f4b38f6d7e5659c2f3478410bf014c01f3
remove(DLBCL_Seurat)

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = DLBCL, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = DLBCL, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = DLBCL, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

#p3 <- MetageneBicorPlot(DLBCL, grouping.var = "conditions", dims.eval = 1:30, 
#                        display.progress = FALSE)
#p3 + geom_smooth(method = 'loess')
DimHeatmap(object = DLBCL, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)

DimHeatmap(object = DLBCL, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

PrintDim(object = DLBCL, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)


<<<<<<< HEAD
#======1.3 QC (skip, ~20k cells were removed)==================================

=======
<<<<<<< HEAD:R/Seurat_setup.R
#======1.3 QC (skip, ~20k cells were removed)==================================
=======
#======1.3 QC ==================================
>>>>>>> 9b8169132a35923308ef95843457f2dbd21713a4:R/Seurat_setup.R
>>>>>>> be7e28f4b38f6d7e5659c2f3478410bf014c01f3
# Run rare non-overlapping filtering
DLBCL <- CalcVarExpRatio(object = DLBCL, reduction.type = "pca",
                               grouping.var = "conditions", dims.use = 1:15)
DLBCL <- SubsetData(DLBCL, subset.name = "var.ratio.pca",accept.low = 0.5)

mito.genes <- grep(pattern = "^mt-", x = rownames(x = DLBCL@data), value = TRUE)
percent.mito <- Matrix::colSums(DLBCL@raw.data[mito.genes, ])/Matrix::colSums(DLBCL@raw.data)
DLBCL <- AddMetaData(object = DLBCL, metadata = percent.mito, col.name = "percent.mito")
DLBCL <- ScaleData(object = DLBCL, genes.use = genes.use, display.progress = FALSE, 
                         vars.to.regress = "percent.mito", do.par = TRUE, num.cores = 4)

#======1.4 align seurat objects (skip, memory overflow)=========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

DLBCL <- AlignSubspace(object = DLBCL, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:20)
#Now we can run a single integrated analysis on all cells!
DLBCL <- RunTSNE(object = DLBCL, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
DLBCL <- FindClusters(object = DLBCL, reduction.type = "cca.aligned", dims.use = 1:20,
                    resolution = 0.8, force.recalc = T, save.SNN = TRUE,
                    n.start = 10, nn.eps = 0.5, print.output = FALSE)



#======1.5 generate seurat objects from Sparse matrix =========================
lnames = load(file = "./data/DLBCL_10.Rda")
lnames
DLBCL <- CreateSeuratObject(raw.data = DLBCL_marged_mat, meta.data = DLBCL_metadata,
                            min.cells = 3,min.genes = 200,names.delim = ".",
                            project = "EC-RW-4311_4262")
DLBCL <- NormalizeData(object = DLBCL, normalization.method = "LogNormalize", scale.factor = 10000)
mito.genes <- grep(pattern = "^MT", x = rownames(x = DLBCL@data), value = TRUE)
percent.mito <- Matrix::colSums(DLBCL@raw.data[mito.genes, ])/Matrix::colSums(DLBCL@raw.data)
DLBCL <- AddMetaData(object = DLBCL, metadata = percent.mito, col.name = "percent.mito")
DLBCL <- ScaleData(object = DLBCL, genes.use = genes.use, display.progress = FALSE, 
                 vars.to.regress = "percent.mito", do.par = TRUE, num.cores = 4)
DLBCL <- RunPCA(object = DLBCL, pc.genes = genes.use, pcs.compute = 100, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = DLBCL, num.pc = 100)
PCHeatmap(DLBCL, pc.use = c(1:3, 70:75), cells.use = 500, do.balanced = TRUE)
DLBCL <- FindClusters(object = DLBCL, reduction.type = "pca", dims.use = 1:75, resolution = 3, 
                    save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)
WGCNA::collectGarbage()
DLBCL <- RunTSNE(object = DLBCL, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", 
               nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", 
               fast_tsne_path = "/Users/yah2014/src/FIt-SNE/bin/fast_tsne", 
               max_iter = 2000)

p1 <- DimPlot(object = DLBCL, reduction.use = "FItSNE", no.legend = F, do.return = TRUE, 
              group.by = "orig.ident",
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = DLBCL, reduction.use = "FItSNE", no.legend = F, group.by = "conditions", 
              do.return = TRUE, vector.friendly = TRUE, pt.size = 0.1) + ggtitle("Tissue") + 
        theme(plot.title = element_text(hjust = 0.5))

plot_grid(p1, p2)
save(DLBCL_marged_mat,DLBCL_metadata,genes.use,DLBCL, file = "./data/DLBCL_10.Rda")


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
save(DLBCL, file = "./data/DLBCL_10_alignment.Rda")
remove(DLBCL.PDX)
remove(DLBCL.primary)
