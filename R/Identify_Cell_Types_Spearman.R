library(Seurat)
library(dplyr)
library(plyr)
library(janitor)
library(pheatmap)
source("./R/Seurat_functions.R")

#====== 4.1 load data  ==========================================
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
                     "12.Unknown")

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

p1 <- TSNEPlot(object = DLBCL, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Primary")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- TSNEPlot(object = DLBCL.PDX, no.legend = TRUE, do.label = TRUE,
               do.return = TRUE, label.size = 6)+
        ggtitle("PDX")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)

#====== 4.2 Organize Immgen data ===========================================
ImmGenDataV1 <- read.csv(file = "./data/ImmGenDataV1.csv")
ImmGenDataV2 <- read.csv(file = "./data/ImmGenDataV2.csv")
ImmGenData <- full_join(ImmGenDataV1, ImmGenDataV2, by = "GeneSymbol")
colnames(ImmGenData) <- sub('X.', '', colnames(ImmGenData)) # remove all X.
ImmGenData <- ImmGenData %>% clean_names()
ImmGenData <- ImmGenData[,-c(1,3,213,214)] # remove probesetid and description

# organize cell type
CellTypes <- colnames(ImmGenData)
CellTypes <- sub('mlp_', 'sc_mlp_', CellTypes) # stem cells
CellTypes <- sub('prob_', 'b_pro_', CellTypes) # B cells
CellTypes <- sub('preb_', 'b_pre_', CellTypes) # B cells
CellTypes <- sub('b1', 'b_b1', CellTypes) # B cells
CellTypes <- sub('pret_', 't_pre_', CellTypes) # T cells
CellTypes <- sub('nkt_', 't_nkt_', CellTypes) # T cells
CellTypes[c(186:198,228:244)] <- paste0("stro_",CellTypes[c(186:198,228:244)]) # stromal cells
CellTypes[248:254] <- paste0("nk_",CellTypes[248:254]) # Innate Lymphocytes
CellTypes <- sub("_",".",CellTypes)
colnames(ImmGenData) <- CellTypes
# reorder
CellTypes <- c(CellTypes[1],sort(CellTypes[-1]))
ImmGenData <- ImmGenData[CellTypes]

# sum up gene expression
ImmGenData$genesymbol <- toupper(ImmGenData$genesymbol)
ImmGenData$genesymbol <- gsub(" ","",ImmGenData$genesymbol)
system.time(ImmGenData <- aggregate(. ~ genesymbol, data=ImmGenData, FUN=sum))
rownames(ImmGenData) <- ImmGenData$genesymbol

# calculate ImmGenData averageExp
ImmGenData_short <- ImmGenData
ImmGenData_short <- ImmGenData_short[,-1]
ImmGenData_short <- as.data.frame(t(ImmGenData_short))
Major_CellTypes <- sub('\\..*', '', rownames(ImmGenData_short))
Full_names <- data.frame("b" = "B_cells",
                          "ba" = "Basophils",
                          "dc" = "Dendritic_cells",
                          "eo" = "Eosinophils",
                          "gn" = "Neutrophils",
                          "mc" = "Mast_cells",
                          "mf" = "Macrophages",
                          "mo" = "Monocytes",
                          "nk" = "Innate_Lymphocytes",
                          "sc" = "Stem_cells",
                          "stro" = "Stromal_cells",
                          "t" = "T_cells",
                          "tgd" = "gd_T_cells")
rownames(Full_names) <- "Major_CellTypes"
ImmGenData_short$Major_CellTypes <- t(Full_names[1,match(Major_CellTypes,names(Full_names))])
ImmGenData_short[,21753:21756]
system.time(ImmGenData_short <- aggregate(. ~ Major_CellTypes, data = ImmGenData_short, FUN=mean))
rownames(ImmGenData_short) <- ImmGenData_short$Major_CellTypes
ImmGenData_short <- as.data.frame(t(ImmGenData_short[,-1]))
ImmGenData_short$genesymbol <- rownames(ImmGenData_short)
head(ImmGenData_short[,(ncol(ImmGenData_short)-3):ncol(ImmGenData_short)])
ImmGenData.summary <- ImmGenData_short
#====== 4.3 Identify Cell Types by Spearman correlation ==================================
Identify_Cell_Types_Spearman <- function(object, gendata,cluster_rows=F,cluster_cols = F,
                                         fontsize_row = 15,fontsize_col = 15,fontsize =20){
        "
        Calculate Average Expression of each ident of seurat object,
        Calculate spearman correlation between above results and provided gendata dataset
        "
        if(class(object) != "seurat") {
                stop(paste("Error : ", object, " is not a seurat object"))
        }
        if(class(gendata) != "data.frame") {
                stop(paste("Error : ", gendata, " is not a data frame"))
        }
        object.AverageExp <- AverageExpression(object)
        object.AverageExp$genesymbol <- toupper(rownames(object.AverageExp))
        object.Exp <- object.AverageExp[,c(ncol(object.AverageExp),
                                       1:(ncol(object.AverageExp)-1))]# move last column to the first
        print("Merge genes expression:")
        table(gendata$genesymbol %in% object.Exp$genesymbol)
        # merge ===============
        object.Exp_gendata <- inner_join(object.Exp, gendata, by = "genesymbol")
        object.Exp_gendata <- object.Exp_gendata[order(object.Exp_gendata$genesymbol),]
        rownames(object.Exp_gendata) <- object.Exp_gendata$genesymbol
        object.Exp_gendata <- object.Exp_gendata[,-1]
        # Spearman correlation primary ==================
        c <- cor(object.Exp_gendata, method="spearman") # or naive_matrix
        diag(c) <-NA
        ident_num <- length(levels(object@ident))
        object_c_gendata <- c[(ident_num+1):nrow(c),1:ident_num]
        pheatmap(object_c_gendata,cex=.9,
                 cluster_rows= cluster_rows,
                 cluster_cols = cluster_cols,
                 fontsize_row = fontsize_row,
                 fontsize_col = fontsize_col,
                 fontsize = fontsize,
                 main = paste("Spearman correlation:",deparse(substitute(object)),
                              "vs.",deparse(substitute(gendata))))
        actual_cell_types <- apply(object_c_gendata, 2, which.max)
        rename_ident <- data.frame("old.ident.ids" = colnames(object_c_gendata),
                                   "new.cluster.ids" = rownames(object_c_gendata)[actual_cell_types])
        print(rename_ident)
        return(rename_ident)
}

# correlate with ImmGenData_short
DLBCL.primary.rename_ident <- Identify_Cell_Types_Spearman(object = DLBCL.primary, gendata = ImmGenData.summary,
                                                    cluster_rows=T,cluster_cols = T)
DLBCL.PDX.rename_ident     <- Identify_Cell_Types_Spearman(object = DLBCL.PDX, gendata = ImmGenData.summary,
                                                    cluster_rows=T,cluster_cols = T)
DLBCL.primary.rename_ident$new.cluster.ids <- c("B cells/Stem cells",
                                                "B cells",
                                                "B cells",
                                                "T/gd T cells",
                                                "B cells",
                                                "B cells",
                                                "B cells/Stem cells",
                                                "T/gd T cells",
                                                "myeloid cells",
                                                "B cells/Stem cells",
                                                "B cells",
                                                "B cells",
                                                "Unknown")

DLBCL.primary@ident <- plyr::mapvalues(x = DLBCL.primary@ident,
                                from = DLBCL.primary.rename_ident$old.ident.ids,
                                to = as.character(DLBCL.primary.rename_ident$new.cluster.ids))

DLBCL.PDX.rename_ident$new.cluster.ids <- c("B cells/Stem cells",
                                            "B cells",
                                            "B cells/Stem cells",
                                            "B cells/Stem cells",
                                            "B cells/Stem cells",
                                            "B cells/Stem cells",
                                            "B cells/Stem cells",
                                            "B cells",
                                            "B cells/Stem cells",
                                            "B cells/Stem cells",
                                            "B cells",
                                            "Unknown",
                                            "B cells/Stem cells")

DLBCL.PDX@ident <- plyr::mapvalues(x = DLBCL.PDX@ident,
                                from = 0:12,
                                to = as.character(DLBCL.PDX.rename_ident$new.cluster.ids))

p1 <- TSNEPlot(object = DLBCL.primary, no.legend = TRUE, do.label = F,
               do.return = TRUE, label.size = 6,colors.use =c("#F8766D","#00BA38",
                                                              "#B79F00","#F564E3","#619CFF"))+
        ggtitle("Primary")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- TSNEPlot(object = DLBCL.PDX, no.legend = TRUE, do.label = F,
               do.return = TRUE, label.size = 6,colors.use =c("#F8766D","#00BA38","#619CFF"))+
        ggtitle("PDX")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)

#===== 4.4  A table with the number of cells of each cluster and subcluster ======
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
table(DLBCL.primary@ident)/DLBCL.primary@data@Dim[2]
table(DLBCL.PDX@ident)/DLBCL.PDX@data@Dim[2]

#===== 4.5  Test Cell Viability using mitochondrial genes expression ======
mt.genes <- grep("^MT-",DLBCL@raw.data@Dimnames[[1]], value = TRUE)
mt.genes <- HumanGenes(DLBCL, mt.genes)
for(i in 1:(length(mt.genes))) print(paste0(i," ",mt.genes[i]))
#mitochondrial genes
p1 <- SingleFeaturePlot.1(object = DLBCL.primary, 
                  feature = "MT-CO3", min.cutoff =  0, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
p2 <- SingleFeaturePlot.1(object = DLBCL.PDX, 
                          feature = "MT-CO3", min.cutoff = 0, 
                  cols.use = c("lightgrey","blue"), pt.size = 0.5)
plot_grid(p1, p2)
