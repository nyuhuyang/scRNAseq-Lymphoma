library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/DLBCL_8.Rda")
lnames

Featureplot <- function(x){
        p <- FeaturePlot(object = DLBCL, 
                         reduction.use = "tsne",
                         features.plot = x, min.cutoff = NA, 
                         cols.use = c("lightgrey","blue"), pt.size = 0.5)
        return(p)
}
#------
Adipocytes <- HumanGenes(DLBCL,c("SLC36A2","P2RX5","MYF5","UCP1","TRIP4","ASCC1"))
Endothelium <- HumanGenes(DLBCL,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                       "Vwf","EMCN","Car4"))
Epithelium <- HumanGenes(DLBCL,c("Epcam","KRT19","KRT5",
                                      "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                                      "Rlbp1","Msln","Upk3b","Lrrn4"))
RPE <- HumanGenes(DLBCL,c("Rpe65","Rlbp1"))
Fibroblast <- HumanGenes(DLBCL,c("FGF1","FGF9","SFRP1"))
#--Hematopoietic----
Hematopoietic <- HumanGenes(DLBCL,c("PTPRC","LAPTM5","SRGN"))
#------Myeloid----
megakaryocytes <-  HumanGenes(DLBCL,c("PPBP","GNG11"))
erythrocyte <-  HumanGenes(DLBCL,c("HBA2","HBB"))
MastCells <- HumanGenes(DLBCL,c("Cma1","Mcpt4","Tpsb2","Cpa3"))
Monocytes <-  HumanGenes(DLBCL,c("LYZ","S100A9","CD14","CCL2","FCGR3A","MS4A7","VMO1"))
Macrophages <- HumanGenes(DLBCL,c("LYZ","CD68","MARCO","Emr1"))
DendriticCells <- HumanGenes(DLBCL,c("Itgax","GPR183","CST3","HLA-DQA1","FCER1A","TSPAN13",
                                     "IL3RA","IGJ"))
Myeloid <-  HumanGenes(DLBCL,c(megakaryocytes,erythrocyte,MastCells,
                               Monocytes,Macrophages,DendriticCells))
#------Lymphoid----
Lymphoid <- HumanGenes(DLBCL,c("Cd19","CD79A","MS4A1",
                                    "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
# T cell
T_Cell <- HumanGenes(DLBCL,c("CD3G","CD3D","CD2","CD8A","IL2RA",
                             "FOXP3","NCAM1","FCGR3A"))
CD4_Naive_T <- HumanGenes(DLBCL,c("CD4","IL7R","GIMAP5","SELL","IL2RG"))
NK <- HumanGenes(DLBCL,c("NKG7","CCL5"))
# B cell
B_StemCell <- HumanGenes(DLBCL,"SPN")
Pre_Pro_B <- HumanGenes(DLBCL,c("CD34","MME","CD38"))
Pro_B <- HumanGenes(DLBCL,c("MME","CD19","SPN","CD38","CD24","IL7","IL3RA"))
Pre_B <- HumanGenes(DLBCL,c("MME","CD19","MS4A1","CD24","CD38","IL7","IL3RA","IL4R"))
Immature_B <- HumanGenes(DLBCL,c("MME","CD19","MS4A1","CR2","CD40","CD24","CD38","IL4R"))
Transitional_B <- HumanGenes(DLBCL,c("CD19","MS4A1","CD5","CR2","CD24","CD38"))
Marginal_zone_B <- HumanGenes(DLBCL,c("CD1C","CD19","MS4A1","CR2","CD27"))
Regulatory_B <- HumanGenes(DLBCL,c("CD1D","CD5","CD19","CR2","CD24"))
Follicular_B <- HumanGenes(DLBCL,c("CD19","MS4A1","CR2","CD22","FCER2","CD24",
                                   "HLA-DRB1","HLA-DQB1","HLA-DRA","HLA-DQA1"))
Activated_B <- HumanGenes(DLBCL,c("CD27","CD19","MS4A1","IL2RA","TNFRSF8","CD69","CD80","CD86","FLT3"))
Germinal_center_B <- HumanGenes(DLBCL,c("MME","CD19","MS4A1","FCER2","CD27","CD38","TNFRSF17"))
Plasma_blast <- HumanGenes(DLBCL,c("CD19","CD38","CD27","TNFRSF17","HLA-DRB1"))
Plasma_cell_long_lived <- HumanGenes(DLBCL,c("CXCR4","CD27","CD38","CD138","CD269"))
Memory_B <- HumanGenes(DLBCL,c("CD19","MS4A1","CD40","CD27","CXCR4","CXCR5","ACKR3"))


Melanocytes <- HumanGenes(DLBCL,c("Pmel","Mlana"))
Mesenchymal <- HumanGenes(DLBCL,c("Pdgfrb","Vim","Has2","Dcn"))
Myelinating_Schwann_cells <- HumanGenes(DLBCL,c("MBP","MPZ"))
Pericytes <- HumanGenes(DLBCL,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                     "Myh11","Mylk","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- HumanGenes(DLBCL,c("Acta2","Myh11"))
Stem_cell <- HumanGenes(DLBCL,c("POU5F1","FUT4","CD34","PROM1","ABCG2","Runx1","ATXN1",
                                     "Nes","NCAM","NGFR"))
Stromal_fibroblasts <- HumanGenes(DLBCL,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- HumanGenes(DLBCL,c("Ihh","Gli1", "Ptch1", "Hhip"))
# Featureplot
Featureplot(Adipocytes) # Adipocytes
Featureplot(Endothelium) # Endothelial Cells
Featureplot(Epithelium[1]) # Epithelium
Featureplot(c(RPE,Melanocytes,Myelinating_Schwann_cells)) # RPE, Melanocytes, Myelinating Schwann cells
Featureplot(Fibroblast) # Fibroblasts
Featureplot(Hematopoietic) # Hematopoietic cells
Featureplot(Myeloid) # Myeloid cells
Featureplot(MastCells)
Featureplot(Monocytes)
Featureplot(Macrophages)
Featureplot(DendriticCells)
Featureplot(Myeloid)
Featureplot(Lymphoid) # Lymphoid cells
# T cell
Featureplot(T_Cell)
Featureplot(CD4_Naive_T)
# B cell
Featureplot(unique(c(B_StemCell,
                     Pre_Pro_B,
                     Pro_B,
                     Pre_B)))
Featureplot(unique(c(Immature_B,
                     Transitional_B,
                     Marginal_zone_B)))
Featureplot(unique(c(Regulatory_B,
                     Activated_B)))
Featureplot(Follicular_B)
Featureplot(Germinal_center_B)
Featureplot(unique(c(Plasma_blast,
                     Plasma_cell_long_lived,
                     Memory_B)))

Featureplot(Mesenchymal) # Mesenchymal cells
Featureplot(Pericytes) # Pericytes
Featureplot(Smooth_muscle_cells)
Featureplot(Stem_cell)
Featureplot(Stromal_fibroblasts)
Featureplot(Neurons)

markers.to.plot <- c()
markers.to.plot <- c(Stem_cell[c(5,6)],Monocytes[1:3],Macrophages,erythrocyte,T_Cell[1:3],NK,
                     CD4_Naive_T[1:3],B_StemCell,Pro_B,Pre_B,Immature_B, Transitional_B,
                     Activated_B,Follicular_B,Germinal_center_B,
                     Plasma_blast,Plasma_cell_long_lived,Memory_B)
markers.to.plot <- c("RUNX1","ATXN1","LYZ","S100A9","CD14","CD68","VMO1","MARCO","EMR1","HBA2",
                     "HBB","CD3G","CD3D","CD2","NKG7","CCL5","CD4","IL7R","GIMAP5",
                     "CD19","CD38","MS4A1","CD27","IL4R",
                     "CD40","CD5","IL2RA","CD69","CD80","CD86","FLT3",
                     "CD22","FCER2","HLA-DRB1","HLA-DQB1","HLA-DRA","HLA-DQA1","TNFRSF17",
                     "CXCR4","CXCR5")
markers.to.plot <- HumanGenes(DLBCL,markers.to.plot, unique =T)
DotPlot(DLBCL, genes.plot = rev(markers.to.plot),
        cols.use = c("cyan","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

DLBCL.subsets <- SplitCells(DLBCL)
DLBCL.PDX <- DLBCL.subsets[[1]]
DLBCL.primary <- DLBCL.subsets[[2]]
DotPlot(DLBCL.primary, genes.plot = rev(markers.to.plot),
        cols.use = c("cyan","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
DotPlot(DLBCL.PDX, genes.plot = rev(markers.to.plot),
        cols.use = c("cyan","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
# Rename ident
table(DLBCL.primary@ident)
idents <- as.data.frame(table(DLBCL.primary@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("B cells 0",
                     "B cells 1",
                     "B cells 2",
                     "B cells 3",
                     "B cells 4",
                     "B cells 5",
                     "B cells 6",
                     "B cells 7",
                     "B cells 8",
                     "B cells 9",
                     "B cells 10",
                     "T & NK cells 11",
                     "B cells 12",
                     "Myeloid cells 13",
                     "T & NK cells 14",
                     "B cells 15")
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
DotPlot(DLBCL.primary, genes.plot = rev(markers.to.plot),
        cols.use = c("cyan","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)
#==============
idents <- as.data.frame(table(DLBCL.PDX@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("B cells 0",
                     "B cells 1",
                     "B cells 2",
                     "B cells 3",
                     "B cells 4",
                     "B cells 5",
                     "B cells 6",
                     "B cells 7",
                     "B cells 8",
                     "T & B cells 9",
                     "T & B cells 10",
                     "T & NK cells 11",
                     "B cells 12",
                     "T & NK cells 13",
                     "T & NK cells 14")
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
                     "T & NK cells",
                     "T & NK cells")
DLBCL.PDX@ident <- plyr::mapvalues(x = DLBCL.PDX@ident,
                                       from = old.ident.ids,
                                       to = new.cluster.ids)
DotPlot(DLBCL.PDX, genes.plot = rev(markers.to.plot),
        cols.use = c("cyan","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

#====== 2.2 Plots ==========================================
p1 <- TSNEPlot(object = DLBCL.primary, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Primary samples")+
        theme(text = element_text(size=15),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- TSNEPlot(object = DLBCL.PDX, no.legend = TRUE, do.label = TRUE,
               do.return = TRUE, label.size = 6)+
        ggtitle("PDX samples")+
        theme(text = element_text(size=15),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p2, p1)
#====== 2.3 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
# Compare clusters for each dataset
SplitTSNEPlot(DLBCL, "conditions")