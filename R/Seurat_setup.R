########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(cowplot)
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
# Load the DLBCL dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here


samples <- c("3119T1", "3119T6","I17_1733", "IP-10927",
             "IP8672", "TR624_M2102", "TR_619_mouse_2080", "TR_619_mouse_2083")
projects <- paste0("EC-RW-",c(rep(4311,4),4262,4262,4311,4311))
conditions <- c("PDX", "PDX","primary","primary",
                "primary","PDX","PDX","PDX")
DLBCL_raw <- list()
DLBCL_Seurat <- list()
for(i in 1:length(samples)){
  DLBCL_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                              samples[i],"/outs/filtered_gene_bc_matrices/hg19/"))
  colnames(DLBCL_raw[[i]]) <- paste0(samples[i],"_",colnames(DLBCL_raw[[i]])) # Must be samples!
  rownames(DLBCL_raw[[i]]) <- sub("hg19_","",rownames(DLBCL_raw[[i]]))
  DLBCL_Seurat[[i]] <- CreateSeuratObject(DLBCL_raw[[i]],min.cells = 3,
                                          min.genes = 200, project = projects[i],
                                          names.delim = ".")
  DLBCL_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
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
g <- lapply(DLBCL_Seurat, function(x) head(rownames(x@hvg.info), 600))
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

#save(DLBCL,genes.use, file = "./data/DLBCL_10_alignment.Rda")
save(DLBCL,genes.use, file = "./data/DLBCL_8_alignment.Rda")
#DLBCL_marged_mat <- as(object = DLBCL@raw.data, "sparseMatrix")
#DLBCL_metadata <- DLBCL@meta.data
#save(DLBCL_marged_mat,DLBCL_metadata,genes.use, file = "./data/DLBCL_8.Rda")
remove(DLBCL_Seurat)

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = DLBCL, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = DLBCL, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = DLBCL, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

p3 <- MetageneBicorPlot(DLBCL, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE) # run on cluster
p3 + geom_smooth(method = 'loess')
DimHeatmap(object = DLBCL, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)

DimHeatmap(object = DLBCL, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

PrintDim(object = DLBCL, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

#======1.3 QC (skip, ~20k cells were removed)==================================

# Run rare non-overlapping filtering
#DLBCL <- CalcVarExpRatio(object = DLBCL, reduction.type = "pca",
#                               grouping.var = "conditions", dims.use = 1:15)
#DLBCL <- SubsetData(DLBCL, subset.name = "var.ratio.pca",accept.low = 0.5)

#mito.genes <- grep(pattern = "^mt-", x = rownames(x = DLBCL@data), value = TRUE)
#percent.mito <- Matrix::colSums(DLBCL@raw.data[mito.genes, ])/Matrix::colSums(DLBCL@raw.data)
#DLBCL <- AddMetaData(object = DLBCL, metadata = percent.mito, col.name = "percent.mito")
#DLBCL <- ScaleData(object = DLBCL, genes.use = genes.use, display.progress = FALSE, 
#                         vars.to.regress = "percent.mito", do.par = TRUE, num.cores = 4)

#======1.4 align seurat objects (memory overflow, run in cluster with 300G memory)=========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
pwd <- "/home/yah2014/Dropbox/Public/Olivier/R/scRNAseq-Lymphoma/"
lnames = load(file = paste0(pwd,"/data/DLBCL_8_alignment.Rda"))
lnames

DLBCL <- AlignSubspace(object = DLBCL, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:20)
#Now we can run a single integrated analysis on all cells!
DLBCL <- RunTSNE(object = DLBCL, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = TRUE)
DLBCL <- FindClusters(object = DLBCL, reduction.type = "cca.aligned", dims.use = 1:20,
                    resolution = 0.8, force.recalc = T, save.SNN = TRUE,
                    n.start = 10, nn.eps = 0.5, print.output = FALSE)
DLBCL <- StashIdent(object = DLBCL, save.name = "ClusterNames_0.8")
DLBCL <- FindVariableGenes(object = DLBCL, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
genes.use <- head(rownames(DLBCL@hvg.info), 1000)

save(DLBCL,genes.use, file = paste0(pwd,"/data/DLBCL_8_alignment.Rda"))
#======1.5 generate seurat objects from Sparse matrix =========================
lnames = load(file = "./data/DLBCL_8_alignment.Rda")
lnames
SplitTSNEPlot(DLBCL)

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
              vector.friendly = F, pt.size = 1) + ggtitle("Sample ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = DLBCL, reduction.use = "FItSNE", no.legend = F, group.by = "conditions", 
              do.return = TRUE, vector.friendly = F, pt.size = 1) + ggtitle("Sample types") + 
        theme(plot.title = element_text(hjust = 0.5)) # vector.friendly
plot_grid(p1, p2)

DLBCL <- StashIdent(object = DLBCL, save.name = "ClusterNames_3")
DLBCL.meta.data <- DLBCL@meta.data
save(DLBCL,genes.use, file = "./data/DLBCL_8.Rda")
save(DLBCL.meta.data,p3,file= "./data/DLBCL_8.metadata.Rda")

#====1.6 correct orig.ident ====
# during CreateSeuratObject, names.delim = "_" incorrectly recorded wrong orig.ident
lnames = load(file = "./data/DLBCL_8.Rda")
lnames
lnames = load(file = "./data/DLBCL_8.metadata.Rda")

DLBCL.meta.data[grepl("TR_619_mouse_2080",rownames(DLBCL.meta.data)),
                "orig.ident"] <- "TR_619_mouse_2080"
DLBCL.meta.data[grepl("TR_619_mouse_2083",rownames(DLBCL.meta.data)),
                "orig.ident"] <- "TR_619_mouse_2083"
DLBCL.meta.data[grepl("I17_1733",rownames(DLBCL.meta.data)),
                "orig.ident"] <- "I17_1733"
table(DLBCL.meta.data$orig.ident)
DLBCL@meta.data$orig.ident <- DLBCL.meta.data$orig.ident
save(DLBCL,genes.use, file = "./data/DLBCL_8.Rda")
p1 <- DimPlot(object = DLBCL, reduction.use = "tsne", no.legend = F, do.return = TRUE, 
              group.by = "orig.ident",
              vector.friendly = F, pt.size = 1) + ggtitle("Sample ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = DLBCL, reduction.use = "tsne", no.legend = F, group.by = "conditions", 
              do.return = TRUE, vector.friendly = F, pt.size = 1) + ggtitle("Sample types") + 
        theme(plot.title = element_text(hjust = 0.5)) # vector.friendly
plot_grid(p1, p2)