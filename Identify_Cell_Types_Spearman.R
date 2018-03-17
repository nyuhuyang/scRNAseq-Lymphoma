library(Seurat)
library(dplyr)
library(plyr)
library(janitor)
library(pheatmap)
source("./R/Seurat_functions.R")

#====== 4.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/DLBCL_alignment.Rda")
lnames
cell.all <- FetchData(DLBCL,"conditions")
cell.primary <- rownames(cell.all)[cell.all$conditions =="primary"]
cell.PDX <- rownames(cell.all)[cell.all$conditions =="PDX"]

DLBCL.primary <- SubsetData(object = DLBCL,
                            cells.use =cell.primary)
DLBCL.PDX <- SubsetData(object = DLBCL,
                        cells.use =cell.PDX)
p1 <- TSNEPlot(object = DLBCL.primary, no.legend = F, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
        ggtitle("Identify conserved cell type markers")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- TSNEPlot(object = DLBCL.PDX, no.legend = F, do.label = TRUE,
               do.return = TRUE, label.size = 6)+
        ggtitle("Identify conserved cell type markers")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)
remove(DLBCL)
#====== 4.2 Organize datasets ==========================================
#=======4.2.1 Organize Immgen data========
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

# calculate ImmGenData averageExp
ImmGenData_short <- ImmGenData
rownames(ImmGenData_short) <- ImmGenData_short$genesymbol
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

#options("experssion" = 5000000)
system.time(ImmGenData_short <- aggregate(. ~ Major_CellTypes, data = ImmGenData_short, FUN=mean))
rownames(ImmGenData_short) <- ImmGenData_short$Major_CellTypes
ImmGenData_short <- as.data.frame(t(ImmGenData_short[,-1]))
ImmGenData_short$genesymbol <- rownames(ImmGenData_short)
head(ImmGenData_short[,(ncol(ImmGenData_short)-3):ncol(ImmGenData_short)])

#=======4.2.1 Calculate AverageExp ============
AverageExp.primary <- AverageExpression(DLBCL.primary)
AverageExp.primary$genesymbol <- toupper(rownames(AverageExp.primary))
AverageExp.primary <- AverageExp.primary[,c(ncol(AverageExp.primary),
                            1:(ncol(AverageExp.primary)-1))]# move last column to the first

AverageExp.PDX <- AverageExpression(DLBCL.PDX)
AverageExp.PDX$genesymbol <- toupper(rownames(AverageExp.PDX))
AverageExp.PDX <- AverageExp.PDX[,c(ncol(AverageExp.PDX),
                                    1:(ncol(AverageExp.PDX)-1))]# move last column to the first
#====== 4.3 Identify cell type ==========================================
# merge
table(ImmGenData_short$genesymbol %in% AverageExp.primary$genesymbol)
table(ImmGenData_short$genesymbol %in% AverageExp.PDX$genesymbol)

DLBCL_ImmGenData.primary <- inner_join(AverageExp.primary, ImmGenData_short, by = "genesymbol")
DLBCL_ImmGenData.PDX <- inner_join(AverageExp.PDX, ImmGenData_short, by = "genesymbol")

DLBCL_ImmGenData.primary <- DLBCL_ImmGenData.primary[order(DLBCL_ImmGenData.primary$genesymbol),]
DLBCL_ImmGenData.PDX <- DLBCL_ImmGenData.PDX[order(DLBCL_ImmGenData.PDX$genesymbol),]

rownames(DLBCL_ImmGenData.primary) <- DLBCL_ImmGenData.primary$genesymbol
rownames(DLBCL_ImmGenData.PDX) <- DLBCL_ImmGenData.PDX$genesymbol

DLBCL_ImmGenData.primary <- DLBCL_ImmGenData.primary[,-1]
DLBCL_ImmGenData.PDX <- DLBCL_ImmGenData.PDX[,-1]

# Spearman correlation primary ==================
c <- cor(DLBCL_ImmGenData.primary, method="spearman") # or naive_matrix
diag(c) <-NA
ident_num <- length(levels(DLBCL.primary@ident))
DLBCL_c_ImmGenData.primary <- c[(ident_num+1):nrow(c),1:ident_num]
pheatmap(DLBCL_c_ImmGenData.primary,cex=.9, #[(ident_num+1):nrow(c),1:ident_num]
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 15,
         fontsize =20,
         main ="Spearman correlation: scRNA-seq vs. ImmGenData correlation")
actual_cell_types <- apply(DLBCL_c_ImmGenData.primary, 2, which.max)
rename_ident <- data.frame("old.ident.ids" = as.numeric(colnames(DLBCL_c_ImmGenData.primary)),
                           "new.cluster.ids" = rownames(DLBCL_c_ImmGenData.primary)[actual_cell_types])
DLBCL.primary@ident <- plyr::mapvalues(x = DLBCL.primary@ident,
                               from = rename_ident$old.ident.ids,
                               to = rename_ident$new.cluster.ids)

DLBCL.primary <- RenameIdent(DLBCL.primary,
                             old.ident.name = rename_ident$old.ident.ids,
                             new.ident.name = rename_ident$new.cluster.ids)
# Spearman correlation primary ==================
c <- cor(DLBCL_ImmGenData.PDX, method="spearman") # or naive_matrix
diag(c) <-NA
ident_num <- length(levels(DLBCL.PDX@ident))
DLBCL_c_ImmGenData.PDX <- c[(ident_num+1):nrow(c),1:ident_num]
pheatmap(DLBCL_c_ImmGenData.PDX,cex=.9, #[(ident_num+1):nrow(c),1:ident_num]
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 15,
         fontsize =20,
         main ="Spearman correlation: scRNA-seq vs. ImmGenData correlation")
actual_cell_types <- apply(DLBCL_c_ImmGenData.PDX, 2, which.max)
rename_ident <- data.frame("old.ident.ids" =colnames(DLBCL_c_ImmGenData.PDX),
                           "new.cluster.ids" = rownames(DLBCL_c_ImmGenData.PDX)[actual_cell_types])
DLBCL.PDX@ident <- plyr::mapvalues(x = DLBCL.PDX@ident,
                                       from = rename_ident$old.ident.ids,
                                       to = rename_ident$new.cluster.ids)

p1 <- TSNEPlot(object = DLBCL.primary, no.legend = F, do.label = TRUE,
               do.return = TRUE, label.size = 6)+
        ggtitle("Identify conserved cell type markers")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- TSNEPlot(object = DLBCL.PDX, no.legend = F, do.label = TRUE,
               do.return = TRUE, label.size = 6)+
        ggtitle("Identify conserved cell type markers")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)