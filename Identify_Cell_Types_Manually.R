library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/DLBCL_alignment.Rda")
lnames

# all marker genes 
marker.genes <- c("Cdh5","Pecam1","Flt1","Vwf","Plvap","Kdr","EMCN","Car4","ptprb",
                  "GPX3","DCN","COL6A1","TIMP3","PDGFRA","Acta2","Myh11","FGF1",
                  "FGF9","SFRP1","PTPRC","LAPTM5","SRGN","Cma1","Mcpt4","Tpsb2",
                  "Cpa3","MS4A7","LYZ","CD14","CCL2","FCGR3A","Emr1","CD68",
                  "MARCO","LYZ","Itgax","FCER1A","CD3G","Cd4","CD62L","IL7R",
                  "IL2RG","SELL","CREM","Cd8a","CREM","SELL","Foxp3","Cd19",
                  "CD79A","MS4A1","CD27","SDC1","IL6R",
                  "SLAMF7","GNLY","Ncr1","KLRD1","NKG7","PPBP","Rbfox3",
                  "Uty","Ddx3y","Xist","Pdgfrb","Cspg4","Anpep","Rgs5",
                  "Myh11","Mylk","Sost","Des","Vtn","Ifitm1","Pdgfrb","Vim",
                  "Has2","Dcn","Eng","CD44","Nt5e","Alcam","Kit","Ly6a","Thy1",
                  "Itgb1","Vcam1","Icam2","CD72","CD2","KRT19","Epcam","KRT5",
                  "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                  "Rlbp1","Msln","Upk3b","Lrrn4","Pmel","Mlana","MBP","MPZ",
                  "SLC36A2","P2RX5","MYF5","UCP1","TRIP4","ASCC1","CD34","PROM1","ABCG2","ATXN1")
marker.genes <- HumanGenes(DLBCL, marker.genes)
for(i in 1:(length(marker.genes))) print(paste0(i," ",marker.genes[i]))
#Endothelial Cells
Endothelial <- marker.genes[1:4] #"FLT1"  "PLVAP" "PTPRB" "GPX3"
FeaturePlot(object = DLBCL, 
            features.plot = Endothelial, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
#fibroblasts
fibroblasts <- marker.genes[c(5:11)] #"COL6A1" "TIMP3"  "PDGFRA" "ACTA2"  "MYH11"  "FGF9"   "SFRP1" 
FeaturePlot(object = DLBCL, 
            features.plot = marker.genes[c(5:11)], min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# Hematopoietic cells
Hematopoietic <- marker.genes[c(12:14)]
myeloid <- marker.genes[c(15:24)]
Monocytes <- marker.genes[c(16:17,19,21:22)]
Tcell <- marker.genes[c(66,25,26:27)]
Bcell <- marker.genes[c(35:41)]
NK <- marker.genes[c(42:46)] #"GNLY"   "NCR1"   "KLRD1"  "NKG7"   "RBFOX3"
FeaturePlot(object = DLBCL, 
            features.plot = Bcell, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
#Pericytes and Mesenchymal Cells
pericytes <- marker.genes[c(48:56)] #"PDGFRB" "ANPEP"  "RGS5"   "MYH11"  "MYLK"   "DES"    "IFITM1" "PDGFRB" "VIM" 
Mesenchymal <- marker.genes[c(57:66)] #"ENG"   "CD44"  "NT5E"  "ALCAM" "THY1"  "ITGB1" "VCAM1" "ICAM2" "CD72"  "CD2"
FeaturePlot(object = DLBCL, 
            features.plot = Mesenchymal, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
#Epithelium
Epithemlium <- marker.genes[c(67:73)] #"EPCAM"   "KRT5"    "MUC1"    "SCGB3A1" "SFTPB"   "FOXJ1"   "UPK3B" 
FeaturePlot(object = DLBCL, 
            features.plot = marker.genes[c(67:73)], min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# RPE ,Melanocytes and Myelinating Schwann cells
RPE <- marker.genes[c(74:79)] #"PMEL"  "MBP"   "MPZ"   "P2RX5" "TRIP4" "ASCC1"
FeaturePlot(object = DLBCL, 
            features.plot = marker.genes[c(74:79)], min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# Hematopoietic Stem cell
stem <- marker.genes[c(80:81)] #"ABCG2" "ATXN1"
FeaturePlot(object = DLBCL, 
            features.plot = stem, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 1-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c("IGJ","CREM","IL2RG","SLAMF7","LYZ","S100A9","CD14",
                     "CD68","MARCO","VMO1","FCGR3A","Itgax","FCER1A","CD3D","CD2",
                     "NKG7","CCL5","Cd4","CD62L","IL7R","GIMAP5","CREM","Cd8a","GNLY",
                     "CCL5","KLRD1","CD34","PROM1","ABCG2","ATXN1","CD79A","MS4A1",
                     "NME1","CACYBP","Cd19","CREM","MIR155HG","HLA-DQA1","HLA-DQA1",
                     "SELL","GPR183","CD27")
markers.to.plot <- HumanGenes(DLBCL,markers.to.plot, unique =T)
sort(markers.to.plot)
sdp <- SplitDotPlotGG(DLBCL, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")

# Rename ident
table(DLBCL@ident)
old.ident.ids <- 0:12
new.cluster.ids <- c("0.B cells",
                     "1.B cells",
                     "2.B cells",
                     "3.T cells",
                     "4.B cells",
                     "5.B cells",
                     "6.B cells",
                     "7.T cells",
                     "8.myeloid cells",
                     "9.B cells",
                     "10.\nB cells",
                     "11.B cells",
                     "12.Unknown")

DLBCL@ident <- plyr::mapvalues(x = DLBCL@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
# DLBCL <- RenameIdentBack(DLBCL)
# How many cells are in each cluster
table(DLBCL@ident)
TSNEPlot(object = DLBCL, no.legend = F, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
  ggtitle("Identify conserved cell type markers")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle

#=====2.2 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
freq_table <- prop.table(x = table(DLBCL@ident, DLBCL@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table


#====== 2.3 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
lnames = load(file = "./data/DLBCL_alignment.Rda")
idents <- as.data.frame(table(DLBCL@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0.B cells",
                     "1.B cells",
                     "2.B cells",
                     "3.T cells",
                     "4.B cells",
                     "5.B cells",
                     "6.B cells",
                     "7.T cells",
                     "8.myeloid cells",
                     "9.B cells",
                     "10.B cells",
                     "11.B cells",
                     "12.T, pDC & plasma cells")

DLBCL@ident <- plyr::mapvalues(x = DLBCL@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
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
