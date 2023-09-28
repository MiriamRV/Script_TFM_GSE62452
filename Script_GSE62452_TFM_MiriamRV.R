# Load libraries
library(GEOquery)
library(affy)
library(tidyverse)
library(biomaRt)
library(limma)
library(openxlsx)
library(plotly)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(factoextra)
library(grid)
library(futile.logger)
library(RColorBrewer)
library(VennDiagram)
library(IMIFA)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(enrichplot)
library(ggupset)
library(gridExtra) 


################################################################################
##                          1. DOWNLOAD THE DATA                              ##
################################################################################
# Download the raw data: .CEL files
getGEOSuppFiles("GSE62452")

# Untar files
untar("./GSE62452/GSE62452_RAW.tar")

# Read raw data
gse_62452 <- ReadAffy(verbose = TRUE)

# Explore data
##gse_62452
#AffyBatch object
#size of arrays=1050x1050 features (78 kb)
#cdf=HuGene-1_0-st-v1 (32321 affyids)
#number of samples=130
#number of genes=32321
#annotation=hugene10stv1
## dim(gse_62452) # Rows 1050 Cols 1050

# Store as array
gse_62452_array <- exprs(gse_62452)
#all(colnames(gse_62452_array) == colnames(gse_62452)) #verifying names
#dim(gse_62452_array) # Rows 1102500 Cols 130


################################################################################
##                          2. PRE-PROCESAMIENTO                              ##
################################################################################

# Colnames follow 2 formats, GSM\\d_Hussain, and GSM\\d
colnames(gse_62452_array) <- gsub("(.+)(_Hussain.+)", "\\1", 
                                  colnames(gse_62452_array)) # Substitute GSM\\d_Hussain
colnames(gse_62452_array) <- gsub("(.+)(_\\d)(.+)", "\\1", 
                                  colnames(gse_62452_array)) # Substitute GSM\\d

# Plot raw data

colores <- c(rep("lightskyblue", 90), rep("royalblue1", 40)) #Color 2 groups

boxplot(gse_62452_array,
        las = 3,  
        ylab = "Expression level", 
        cex.axis = 0.6, 
        outline = FALSE,
        col = colores)

legend("topleft", legend = c("Dataset1", "Dataset2"),
       fill = c("lightskyblue", "royalblue1"))

# a) Correct Background, Normalize intra-array & transform from probe to probeset
corregido_62452 <- expresso(gse_62452, 
                            bg.correct = TRUE,
                            bgcorrect.method = 'rma',
                            normalize = TRUE, 
                            normalize.method = NULL, 
                            pmcorrect.method = 'pmonly', 
                            summary.method = 'avgdiff')

expression_level <- exprs(corregido_62452)
colnames(expression_level) <- colnames(gse_62452_array)
#class(expression_level) #"matrix" "array"
#dim(expression_level) # Rows 32321 (probes) cols 130 (samples)

# Plot background
boxplot(expression_level, 
        las = 3,  
        ylab = "Fluorescencia", 
        cex.axis = 0.6, 
        outline = FALSE,
        col = colores)

legend("topleft", legend = c("Dataset1", "Dataset2"),
       fill = c("lightskyblue", "royalblue1"))

# b) From probeset to gen: ANNOTATION 
# Download annotation
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
# Check genome version in ensembl
## searchDatasets(mart = mart, pattern = "hsapiens")  #version GRCh38.p13
annot <- getBM(attributes = c("affy_hugene_1_0_st_v1", "hgnc_symbol"),
               filters = "affy_hugene_1_0_st_v1",
               values = rownames(expression_level),
               mart = mart,
               uniqueRows=TRUE)

colnames(annot) <- c('sondas', 'genes')

#dim(annot) # Rows 36687 cols 2
#class(annot) # dataframe

# Remove non-associated-probesets 
# Sustitute space for NA
df_NA <- annot[annot == ''] <- NA #Add NA to annot
#sum(is.na(annot)) #3301

# Remove NA
df <- na.omit(annot)
colnames(df) <- c('sondas', 'genes')
#dim(df) # Rows 33386 cols 2
#sum(is.na(df)) # 0

# Ckeck duplicates
sum(duplicated(df)) # 0

# resumo con la mediana los valores de expresión de las probsetes asociadas al mismo gen:
mat_annot = do.call("rbind", lapply(unique(df$genes), function(x){
  
  sondas <- as.character(df$sondas[df$genes == x])
  
  if (x!="’"){summ_expr <- apply(matrix(expression_level[sondas,], 
                                        ncol = ncol(expression_level)), 2, median)}
  
}))

rownames(mat_annot) <- unique(df$genes)
colnames(mat_annot) <- colnames(gse_62452_array)


sum(is.na(mat_annot)) # Number of probesets not matching any gene = 0

# log2 transformation
mat_log <- log2(mat_annot)

# c) Normalization between arrays
mat_norm <- normalizeBetweenArrays(mat_log, method = "quantile")
rownames(mat_norm) <- rownames(mat_log)
colnames(mat_norm) <- colnames(mat_log)

#dim(mat_annot) # Rows 25380 (genes)  cols 130 (samples)

boxplot(mat_norm, 
        las = 3, 
        ylab = "Expression level", 
        cex.axis = 0.6, 
        outline = FALSE,
        col = colores)

legend("topleft", legend = c("Dataset1", "Dataset2"),
       fill = c("lightskyblue", "royalblue1"))

##################### PRINCIPAL COMPONENT ANALYSIS (PCA) #######################

# Obtain metadata
GSE <- getGEO("GSE62452")
Metadata <- pData(phenoData(GSE[[1]]))
Metadata <- select(Metadata, 'title', 'geo_accession','grading:ch1', 
                   'source_name_ch1', 'Stage:ch1','survival months:ch1',
                   'survival status:ch1','tissue:ch1') # Select the most importants



## Modify Metadata (add 2 new columns to facilitate the analysis) 
### 1. Grouping the stages 
factor(Metadata$`Stage:ch1`)
Metadata$stage_group <- factor(Metadata$`Stage:ch1`, 
                               labels = c("III", "I", "II", "II", "III", 
                                          rep("IV", 3)))
Metadata$stage_group <- factor(Metadata$stage_group, ordered = TRUE, 
                               levels = c("I", "II", "III", 'IV'))
#levels(Metadata$stage_group) #"I" "II" "III" "IV" 

### 2. Add a new column with the patients (45 patients/Paired samples)
Metadata$Patients = NA
Metadata$Patients[1:90] <- rep(c(1:45), each = 2) 

### Change colnames metadata
#names(Metadata)
colnames(Metadata)[1] <- 'Sample'
colnames(Metadata)[2] <- 'geo_accession'
colnames(Metadata)[3] <- 'grading'
colnames(Metadata)[4] <- 'Source'
colnames(Metadata)[5] <- 'stage'
colnames(Metadata)[6] <- 'Survival_status'
colnames(Metadata)[7] <- 'Survival_status2'
colnames(Metadata)[8] <- 'Tissue'
row.names(Metadata) = NULL 

### Replace "?" by "NA" in Metadata$Survival_status and Survival_status2
Metadata$Survival_status[Metadata$Survival_status == "?"] <- NA
Metadata$Survival_status2[Metadata$Survival_status2 == "?"] <- NA

### Convert Metadata$Survival_status in class numeric
Metadata$Survival_status <- as.numeric(Metadata$Survival_status)

## Create an excel file with the Medatada
packageDescription("openxlsx")$Version
BiocManager::install("openxlsx")

write.xlsx(Metadata, "metadataGSE62452.xlsx")

## Convert the normalised matrix to dataframe
matrix_expression <- as.data.frame(mat_norm)

## a) PCA ALL: tumor tissue and adjacent pancreatic non-tumor (130 samples)
res_pcaALL <- prcomp(t(matrix_expression), scale. = TRUE)
summary(res_pcaALL)

plot_pca12_ALL <- autoplot(res_pcaALL, x = 1, y = 2, 
                           data = Metadata, 
                           colour = 'Tissue',
                           shape = 'stage_group',
                           size = 3,
                           frame = TRUE,
                           cex.lab = 6,
                           label.size = 8,
                           frame.type = "norm") +
                           theme_bw()
ggplotly(plot_pca12_ALL)
ggsave("PCA_TissueALL12-2.jpg")

plot_pca13_ALL <- autoplot(res_pcaALL, x = 1, y = 3, 
                           data = Metadata, 
                           colour = 'Tissue',
                           shape = 'stage_group',
                           size = 3,
                           frame = TRUE,
                           frame.type = "norm")  +
                           theme_bw()
ggplotly(plot_pca13_ALL)
ggsave("PCA_TissueALL13.jpg")

plot_pca23_ALL <- autoplot(res_pcaALL, x = 2, y = 3, 
                           data = Metadata, 
                           colour = 'Tissue',
                           shape = 'stage_group',
                           size = 3, 
                           frame = TRUE,
                           frame.type = "norm") +
                           theme_bw()
ggplotly(plot_pca23_ALL)
ggsave("PCA_TissueALL23.jpg")

fviz_eig(res_pcaALL, addlabels = TRUE, ## Scree plot ALL
         linecolor = "Red", ylim = c(0, 50))

## b) PAIRED SAMPLES: 45 tumors & 45 adjacent pancreatic non-tumor (90 samples)
res_pca_Paired <- prcomp(t(matrix_expression[,1:90]), scale. = TRUE)
summary(res_pca_Paired)


plot_pca12_Paired <- autoplot(res_pca_Paired, x = 1, y = 2, 
                              data = Metadata[1:90,], 
                              colour = 'Tissue',
                              size = 3,
                              frame = TRUE,
                              frame.type = "norm",
                              label.size = 4, 
                              label = TRUE,
                              label.label = "Patients",
                              label.colour = 'black') +
                              theme_bw()

ggplotly(plot_pca12_Paired)
ggsave("plot_pca12_Paired.jpg")

plot_pca13_Paired <- autoplot(res_pca_Paired, 
                              x = 1, y = 3, 
                              data = Metadata[1:90,], 
                              colour = 'Tissue',
                              size = 3,
                              frame = TRUE,
                              frame.type = 'norm',
                              label = TRUE,
                              label.size = 4,
                              label.label = "Patients",
                              label.colour = 'black') + 
                              theme_bw()
ggplotly(plot_pca13_Paired)
ggsave("plot_pca13_Paired.jpg")

plot_pca23_Paired <- autoplot(res_pca_Paired, 
                              x = 2, y = 3, 
                              data = Metadata[1:90,], 
                              colour = 'Tissue',
                              size = 3,
                              frame = TRUE,
                              frame.type = 'norm',
                              label.size = 4, 
                              label = TRUE,
                              label.label = "Patients",
                              label.colour = 'black') +
                              theme_bw()
ggplotly(plot_pca23_Paired)
ggsave("plot_pca23_Paired.jpg")

fviz_eig(res_pca_Paired, addlabels = TRUE,  ## Scree plot Paired
         linecolor = "Red", ylim = c(0, 50))

# c) PCA TUMORS (69 samples)
res_pca <- prcomp(t(matrix_expression[, colnames(matrix_expression) %in%
                                        Metadata[Metadata$Tissue == "Pancreatic tumor", ]$geo_accession]))
summary(res_pca)

## PCA according to stages (ONLY TUMORS!)

# Define colour pallete
palette = c('#990F02', '#005A92', '#00A86B', '#FDEE00')

plot_pca12_stagesT <- autoplot(res_pca, x = 1, y = 2, data = Metadata %>% 
                                 filter(`Tissue` == 'Pancreatic tumor') , 
                               colour = 'stage_group',
                               size = 2,
                               frame = TRUE) +
                               theme_bw() +
                               scale_colour_manual(values = palette)

ggplotly(plot_pca12_stagesT)
ggsave("PCA_stagesT12.jpg")

plot_pca23_stagesT <- autoplot(res_pca, x = 2, y = 3, data = Metadata %>% 
                                 filter(`Tissue` == 'Pancreatic tumor'), 
                               colour = 'stage_group',
                               size = 2,
                               frame = TRUE) +
                               theme_bw() +
                               scale_colour_manual(values = palette)
                               
ggplotly(plot_pca23_stagesT)
ggsave("PCA_stagesT23.jpg")

plot_pca13_stagesT <- autoplot(res_pca, x = 1, y = 3, data = Metadata %>% 
                                 filter(`Tissue` == 'Pancreatic tumor'), 
                               colour = 'stage_group',
                               size = 2,
                               frame = TRUE) +
                               theme_bw() +
                               scale_colour_manual(values = palette)
ggplotly(plot_pca13_stagesT)
ggsave("PCA_stagesT13.jpg")

fviz_eig(res_pca, addlabels = TRUE, ## Scree plot Tumors
         linecolor = "Red", ylim = c(0, 50))

## PCA according to the grading (ONLY TUMORS!)
plot_pca12_gradingT <- autoplot(res_pca, x = 1, y = 2, data = Metadata %>%
                                  filter(`Tissue` == 'Pancreatic tumor'), 
                                colour = 'grading',
                                size = 3) + theme_bw()
ggplotly(plot_pca12_gradingT)
ggsave("PCA_gradingT12.jpg")

plot_pca13_gradingT <- autoplot(res_pca, x = 1, y = 3, data = Metadata %>% 
                                  filter(`Tissue` == 'Pancreatic tumor'), 
                                colour = 'grading',
                                size = 3) + theme_bw()
ggplotly(plot_pca13_gradingT)
ggsave("PCA_gradingT13.jpg")

plot_pca23_gradingT <- autoplot(res_pca, x = 2, y = 3, data = Metadata %>% 
                                  filter(`Tissue` == 'Pancreatic tumor'), 
                                colour = 'grading',
                                size = 3) + theme_bw()
ggplotly(plot_pca23_gradingT)
ggsave("PCA_gradingT23.jpg")

## PCA according to survival (ONLY TUMORS!)
plot_pca12_survivalT <- autoplot(res_pca, x = 1, y = 2, data = Metadata %>% 
                                   filter(`Tissue` == 'Pancreatic tumor'), 
                                 colour = 'Survival_status',
                                 shape = 'stage_group',
                                 size = 3) +
                                 scale_colour_gradient(low = 'blue', 
                                                       high = 'red') +
                                 theme_bw()

ggplotly(plot_pca12_survivalT)
ggsave("PCA_survivalT.jpg")

plot_pca13_survivalT <- autoplot(res_pca, x = 1, y = 3, data = Metadata %>% 
                                   filter(`Tissue` == 'Pancreatic tumor'), 
                                 colour = 'Survival_status',
                                 shape = 'stage_group',
                                 size = 3) +
                                 scale_colour_gradient(low = 'blue', 
                                                       high = 'red') +
                                 theme_bw()

ggplotly(plot_pca13_survivalT)
ggsave("PCA_survivalT13.jpg")

plot_pca23_survivalT <- autoplot(res_pca, x = 2, y = 3, data = Metadata %>% 
                                   filter(`Tissue` == 'Pancreatic tumor'), 
                                 colour = 'Survival_status',
                                 shape = 'stage_group',
                                 size = 3) +                                 
                                 scale_colour_gradient(low = 'blue', 
                                                       high = 'red') +
                                 theme_bw()

ggplotly(plot_pca23_survivalT)
ggsave("PCA_survivalT23.jpg")

################################################################################
##                    3. DIFFERENTIAL GENE EXPRESSION (DGE)                   ##
################################################################################

# 3.1. Analysis tumor tissue and adjacent pancreatic non-tumor (130 samples)

## 3.1.1. ALL samples: Analysis tumor and adjacent tissue (130 samples)

### Design
TvC_design <- model.matrix(~0 + Metadata$Tissue, data = Metadata)
colnames(TvC_design) <- c("Tumor", "Control")
rownames(TvC_design) <- colnames(matrix_expression)

### Lineal model ALL
fit_TvC <- lmFit(matrix_expression, TvC_design)

### Compare Tumor Vs control tissue
contrast.matrix_TvC <- makeContrasts(Tumor-Control, levels = TvC_design)
fit_TvC2 <- contrasts.fit(fit_TvC, contrast.matrix_TvC)
fit_TvC2_results <- eBayes(fit_TvC2)
Results_designTvC <- decideTests(fit_TvC2_results)
summary(decideTests(fit_TvC2_results))  #Tumor-Control: Down 3584 UP 3558

### A list of top genes differential expressed 
TvC <- topTable(fit_TvC2_results, number = 25344, coef = 1, adjust = "BH")

###
TvC$diffexpressed = "No sig."

TvC$diffexpressed[TvC$logFC > 1.5 & TvC$adj.P.Val < 0.05] <- "UP"
TvC$diffexpressed[TvC$logFC < -1.5 & TvC$adj.P.Val < 0.05] <- "DOWN"

sum(with(TvC,diffexpressed == "UP")) # 62 genes
sum(with(TvC,diffexpressed == "DOWN")) # 67 genes

TvC$genes <- row.names.data.frame(TvC)

TvC$delabel <- NA
TvC$delabel[TvC$diffexpressed != "No sig."] <- TvC$genes[TvC$diffexpressed != "No sig."]


top <- 10 #Top 10 DEG
top_genes_TvC <- bind_rows(
  TvC %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top),
  TvC %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top))
colnames(top_genes_TvC)

Table_TvC <- top_genes_TvC[,c(8,1,4,5,7)]
write.xlsx(Table_TvC, "table_TvC.xlsx")

DEG <- 25344 #All DEG
DEG_genes_TvC <- bind_rows(
  TvC %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG),
  TvC %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG))

### Volcano plot with ggplot2 ALL samples
ggplot(data = TvC, aes(x = logFC, y = -log10(adj.P.Val), 
                       col = diffexpressed, label = delabel)) +
  geom_point(size = 2) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "darkgrey", "red")) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col ="black", linetype = "dashed") +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 20)) 
    
ggsave("volcano_ALL-2.pdf")

### Heatmap
hmcol = greenred(1000)
grouped <- arrange(Metadata, Tissue) #Group by Tissue
conds = grouped$Tissue
conds = as.factor(conds)
cond_colours = rainbow(length(unique(conds)), 
                       alpha = 1)[as.factor(conds)] 
sig_genes_ALL<- subset(matrix_expression, 
                          subset = rownames(matrix_expression) %in% TvC$delabel)
col_order_ALL <- c(grouped$geo_accession)
sig_genes_ALL<- sig_genes_ALL[ ,col_order_ALL]
inter = unique(grouped$Tissue)
a = factor(conds, levels = inter) 

pdf("heatmap_ALL_Group.pdf") ## Tumors vs Control samples
heatmap(as.matrix(sig_genes_ALL), 
        col = hmcol, 
        cexRow = 0.5, 
        cexCol = 0.1,
        Colv = NA,
        ColSideColors = cond_colours)
legend("topleft", 
       levels(a),
       fill = unique(cond_colours), 
       cex = 0.5)
heat_legend(as.matrix(sig_genes_ALL), 
            hmcol,
            breaks = c(-2, 0, 2),
            cex.lab = 0.3)
dev.off()

pdf("heatmap_ALL_Cluster.pdf") ## Clustering 
heatmap(as.matrix(sig_genes_ALL), 
        col = hmcol, 
        cexRow = 0.5, 
        cexCol = 0.1,
        ColSideColors = cond_colours)
heat_legend(as.matrix(sig_genes_ALL), 
            hmcol,
            cex.lab = 0.3,
            size = 3)
legend("topleft",
       levels(a), 
       fill = unique(cond_colours), 
       cex = 0.3)
dev.off()

## 3.2.2. PAIRED SAMPLES: 45 tumors & 45 adjacent pancreatic non-tumor (90 samples)

### Design paired
Patients <- factor(Metadata[1:90,]$Patients)
Tissue <- factor(Metadata[1:90,]$Tissue, 
                 levels = c("Pancreatic tumor", "adjacent pancreatic non-tumor"))
paired_design <- model.matrix(~Patients+Tissue) 
rownames(paired_design) <- colnames(matrix_expression[,1:90])

### Lineal Model paired
fit_paired <- lmFit(matrix_expression[,1:90], paired_design)
fit_paired2 <- eBayes(fit_paired)

Resultados_designPaired <- decideTests(fit_paired2)
summary(decideTests(fit_paired2))

### A list of top genes differential expressed 
Paired_results <- topTable(fit_paired2,  number = 25344, coef = 46, adjust = "BH")

###
Paired_results$diffexpressed = "No sig."

Paired_results$diffexpressed[Paired_results$logFC > 1.5 & Paired_results$adj.P.Val < 0.05] <- "UP"
Paired_results$diffexpressed[Paired_results$logFC < -1.5 & Paired_results$adj.P.Val < 0.05] <- "DOWN"

sum(with(Paired_results,diffexpressed== "UP")) # 91 genes
sum(with(Paired_results,diffexpressed== "DOWN")) # 92 genes

Paired_results$genes <- row.names.data.frame(Paired_results)

Paired_results$delabel <- NA
Paired_results$delabel[Paired_results$diffexpressed != "No sig."] <- Paired_results$genes[Paired_results$diffexpressed != "No sig."]

top_genes_Paired <- bind_rows(  #Top10 DEG
  Paired_results %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top),
  Paired_results %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top))

Table_Paired <- top_genes_Paired[,c(8,1,4,5,7)]
write.xlsx(Table_Paired, "table_Paired.xlsx")


DEG_genes_Paired <- bind_rows(   #All DEG
  Paired_results %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG),
  Paired_results %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG))

### Volcano plot with ggplot2 Paired samples
ggplot(data = Paired_results, aes(x = logFC, y = -log10(adj.P.Val), 
                                  col = diffexpressed, label = delabel)) +
  geom_point(size = 2) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "darkgrey", "red")) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col ="black", linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12)) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 15))+
  theme(legend.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) 

ggsave("volcano_paired.pdf")

### Heatmap Tumor Vs Control
grouped_Paired <- Metadata[1:90,]
grouped_Paired <- arrange(grouped_Paired, Tissue)
conds_Paired = grouped_Paired$Tissue
conds_Paired = as.factor(conds_Paired)
sig_genes_Paired <- subset(matrix_expression[,1:90], 
                           subset = rownames(matrix_expression) %in% Paired_results$delabel)
col_order_Paired <- c(grouped_Paired$geo_accession)
sig_genes_Paired <- sig_genes_Paired[ ,col_order_Paired]
cond_colours_Paired = rainbow(length(unique(conds_Paired)), 
                              alpha = 1)[as.factor(conds_Paired)] 
inter_Paired = unique(grouped_Paired$Tissue)
b = factor(conds_Paired, levels = inter_Paired)

pdf("heatmap_Paired_Group.pdf")
heatmap(as.matrix(sig_genes_Paired), 
        col = hmcol, 
        cexRow = 0.5, 
        cexCol = 0.5,
        Colv = NA,
        ColSideColors = cond_colours_Paired)
heat_legend(as.matrix(sig_genes_Paired), 
            hmcol,
            breaks = c(-2, 0, 2),
            cex.lab = 0.4)
legend("topleft", 
       levels(b), 
       fill = unique(cond_colours_Paired), 
       cex = 0.3)
dev.off()

pdf("heatmap_Paired_Cluster.pdf")
heatmap(as.matrix(sig_genes_Paired), 
        col = hmcol, 
        cexRow = 0.5, 
        cexCol = 0.5,
        ColSideColors = cond_colours_Paired)
heat_legend(as.matrix(sig_genes_Paired), 
            hmcol, 
            cex.lab = 0.4)
legend("topleft", 
       levels(b), 
       fill = unique(cond_colours_Paired), 
       cex = 0.3)
dev.off()

### Heatmap Patients
conds_Patients = Metadata$Patients[1:90]
conds_Patients = as.factor(conds_Patients)
cond_colours_Patients = rainbow(length(unique(conds_Patients)), 
                              alpha = 1)[as.factor(conds_Patients)] 
sig_genes_Paired_Patients <- subset(matrix_expression[,1:90], 
                           subset = rownames(matrix_expression) %in% Paired_results$delabel)
inter_Patients = unique(Metadata$Patients[1:90])
b2 = factor(conds_Patients, levels = inter_Patients) 

pdf("heatmap_Patients.pdf")
heatmap(as.matrix(sig_genes_Paired_Patients), 
        col = hmcol, 
        cexRow = 0.5, 
        cexCol = 0.5,
        Colv = NA,
        ColSideColors = cond_colours_Patients)
heat_legend(as.matrix(sig_genes_Paired_Patients), 
            hmcol, 
            cex.lab = 0.4)
legend("left",
       yjust = 0,
       xjust = 0,
       levels(b2), 
       fill = unique(cond_colours_Patients), 
       cex = 0.3)
dev.off()

pdf("heatmap_Patients_exp.pdf")
heatmap(as.matrix(sig_genes_Paired_Patients), 
        col = hmcol, 
        cexRow = 0.5, 
        cexCol = 0.5,
        ColSideColors = cond_colours_Patients, cond_colours_Paired)
dev.off()

### Venn Diagram Tumors vs Control (ALL & Paired)
deg2 = c(TvC$genes)
deg3 = c(Paired_results$genes)

venn.diagram(
  x = list(deg2, deg3),
  category.names = c("ALL", "Paired"),
  filename = "Venndiagram_TvC.png",
  output = TRUE,
  lwd = 2,  # Circles
  lty = 'blank',
  fill = c("lightgreen", "lightsalmon"),
  col = c("lightgreen", "lightsalmon"),
  cex = 0.6,  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("green", "salmon"),
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055))


# 3.2. TUMORS

## DESIGN
tumors_design <- model.matrix(~0 + (Metadata %>% filter(`Tissue` == 'Pancreatic tumor'))$stage_group)
colnames(tumors_design) <- c("I", "II", "III", "IV")
rownames(tumors_design) <- rownames((Metadata %>% filter(`Tissue` == 'Pancreatic tumor'))$geo_accession)

## Create lineal model
fit <- lmFit(matrix_expression[, colnames(matrix_expression) %in% 
                                  Metadata[Metadata$Tissue == "Pancreatic tumor", ]$geo_accession], tumors_design)

## Compare by pairs the different stages
contrast.matrix <- makeContrasts(II-I, III-I, IV-I, III-II, IV-II, IV-III, levels = tumors_design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit_results <- eBayes(fit2)
Results_design <- decideTests(fit_results)
summary(decideTests(fit_results))

## Compare by pairs only the interesting stages
contrast.matrix2 <- makeContrasts(II-I, III-I, levels = tumors_design)
fit3 <- contrasts.fit(fit, contrast.matrix2)
fit_results_new <- eBayes(fit3)
Results_design2 <- decideTests(fit_results_new)
summary(decideTests(fit_results_new))

## List of top genes differential expressed in group2 versus group1
II_I <- topTable(fit_results, number = 25344, coef = 1, adjust.method = "BH")
III_I <- topTable(fit_results, number = 25344, coef = 2, adjust.method = "BH")
IV_I <- topTable(fit_results, number = 25344, coef = 3, adjust.method = "BH")
III_II <- topTable(fit_results, number = 25344, coef = 4, adjust.method = "BH")
IV_II <- topTable(fit_results, number = 25344, coef = 5, adjust.method = "BH")
IV_III <- topTable(fit_results, number = 25344, coef = 6, adjust.method = "BH")

## Comparation II_I
II_I$diffexpressed = "No sig."

II_I$diffexpressed[II_I$logFC > 0 & II_I$adj.P.Val < 0.05] <- "UP"
II_I$diffexpressed[II_I$logFC < 0 & II_I$adj.P.Val < 0.05] <- "DOWN"


sum(with(II_I,diffexpressed == "UP")) # 1 genes
sum(with(II_I,diffexpressed == "DOWN")) # 3 genes

II_I$genes <- row.names.data.frame(II_I)

II_I$delabel <- NA
II_I$delabel[II_I$diffexpressed != "No sig."] <- II_I$genes[II_I$diffexpressed != "No sig."]


DEG_genes_II_I <- bind_rows( #All DEG
  II_I %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG),
  II_I %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG))

Table_Tumors_II_I <- top_genes_II_I[,c(8,1,4,5,7)]
write.xlsx(Table_Tumors_II_I, "table_Tumors.xlsx")

## Volcano plot with ggplot2 II_I
ggplot(data = II_I, aes(x = logFC, y = -log10(adj.P.Val), 
                        col = diffexpressed, label = delabel)) +
  geom_point(size = 2) +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "darkgrey", "red")) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 15))+
  theme(legend.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) 

ggsave("volcano_II_I.pdf")

## Comparation III_I
III_I$diffexpressed = "No sig."

III_I$diffexpressed[III_I$logFC > 1 & III_I$adj.P.Val < 0.05] <- "UP"
III_I$diffexpressed[III_I$logFC < -1 & III_I$adj.P.Val < 0.05] <- "DOWN"


sum(with(III_I,diffexpressed == "UP")) # 18 genes
sum(with(III_I,diffexpressed == "DOWN")) # 7 genes

III_I$genes <- row.names.data.frame(III_I)

III_I$delabel <- NA
III_I$delabel[III_I$diffexpressed != "No sig."] <- III_I$genes[III_I$diffexpressed != "No sig."]

## Volcano plot with ggplot2 III_I
ggplot(data = III_I, aes(x = logFC, y = -log10(adj.P.Val), 
                         col = diffexpressed, label = delabel)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "darkgrey", "red")) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col ="black", linetype = "dashed") +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 15))+
  theme(legend.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) 

ggsave("volcano_III_I.pdf")

top_genes_III_I <- bind_rows(  #Top 10 DEG
  III_I %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top),
  III_I %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top))

Table_Tumors <- top_genes_III_I[,c(8,1,4,5,7)]
write.xlsx(Table_Tumors, "table_Tumors.xlsx")

DEG_genes_III_I <- bind_rows( #All DEG
  III_I %>% 
    filter(diffexpressed == 'UP') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG),
  III_I %>% 
    filter(diffexpressed == 'DOWN') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(DEG))

# Venn Diagram

# a) All genes
vennDiagram(Results_design2, 
            circle.col = c("blue", "green"))

# b ) DEG
DEG0 = c(DEG_genes_II_I$genes)
DEG1 = c(DEG_genes_III_I$genes)

venn.diagram(
  x = list(DEG0, DEG1), 
  category.names = c("II_I", "III_I"),
  filename = "Venndiagram_DEG_Tumors3.png",
  output = TRUE,
  lwd = 2,  # Circles
  lty = 'blank',
  fill = c("blue", "green"),
  alpha = 0.3,
  cex = 0.6,  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 3,  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer", 
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("darkblue", "darkgreen"))
 ??venn.diagram

## Comparation IV_I
IV_I$diffexpressed = "No sig."

IV_I$diffexpressed[IV_I$logFC > 0 & IV_I$adj.P.Val < 0.05] <- "UP"
IV_I$diffexpressed[IV_I$logFC < 0 & IV_I$adj.P.Val < 0.05] <- "DOWN"


sum(with(IV_I,diffexpressed == "UP")) # 0 genes
sum(with(IV_I,diffexpressed == "DOWN")) # 0 genes

IV_I$genes <- row.names.data.frame(IV_I)

IV_I$delabel <- NA
IV_I$delabel[IV_I$diffexpressed != "No sig."] <- IV_I$genes[IV_I$diffexpressed != "No sig."]

## Comparation III_II
III_II$diffexpressed = "No sig."

III_II$diffexpressed[III_II$logFC > 0 & III_II$adj.P.Val < 0.05] <- "UP"
III_II$diffexpressed[III_II$logFC < 0 & III_II$adj.P.Val < 0.05] <- "DOWN"


sum(with(III_II,diffexpressed == "UP")) # 0 genes
sum(with(III_II,diffexpressed == "DOWN")) # 0 genes

III_II$genes <- row.names.data.frame(III_II)

III_II$delabel <- NA
III_II$delabel[III_II$diffexpressed != "No sig."] <- III_II$genes[III_II$diffexpressed != "No sig."]

## Comparation IV_II
IV_II$diffexpressed = "No sig."

IV_II$diffexpressed[IV_II$logFC > 0 & IV_II$adj.P.Val < 0.05] <- "UP"
IV_II$diffexpressed[IV_II$logFC < 0 & IV_II$adj.P.Val < 0.05] <- "DOWN"


sum(with(IV_II,diffexpressed == "UP")) # 0 genes
sum(with(IV_II,diffexpressed == "DOWN")) # 0 genes

IV_II$genes <- row.names.data.frame(IV_II)

IV_II$delabel <- NA
IV_II$delabel[IV_II$diffexpressed != "No sig."] <- IV_II$genes[IV_II$diffexpressed != "No sig."]

## Comparation IV_III
IV_III$diffexpressed = "No sig."

IV_III$diffexpressed[IV_III$logFC > 0 & IV_III$adj.P.Val < 0.05] <- "UP"
IV_III$diffexpressed[IV_III$logFC < 0 & IV_III$adj.P.Val < 0.05] <- "DOWN"


sum(with(IV_III,diffexpressed == "UP")) # 0 genes
sum(with(IV_III,diffexpressed == "DOWN")) # 0 genes

IV_III$genes <- row.names.data.frame(IV_III)

IV_III$delabel <- NA
IV_III$delabel[IV_III$diffexpressed != "No sig."] <- IV_III$genes[IV_III$diffexpressed != "No sig."]

# Venn diagram for all the conditions: ALL, Paired and Tumors 

# a) All genes 
deg1 = c(III_I$genes)
deg2 = c(TvC$genes)
deg3 = c(Paired_results$genes)
myCol <- brewer.pal(3, "Pastel2") # Prepare a palette of 3 colors

venn.diagram(
  x = list(deg1, deg2, deg3), 
  category.names = c("Tumores", "Todos", "Pareados"),
  filename = "Venndiagram_31.png",
  output = TRUE,
          lwd = 2,  # Circles
          lty = 'blank',
          fill = myCol,
          cex = 0.6,  # Numbers
          fontface = "bold",
          fontfamily = "sans",
          cat.cex = 2,  # Set names
          cat.fontface = "bold",
          cat.default.pos = "outer",
          cat.col = c("darkgreen", "darkred", "darkblue"),
          cat.pos = c(-27, 27, 135),
          cat.dist = c(0.055, 0.055, 0.085),
          cat.fontfamily = "sans",
          rotation = 1)

# b) DEG
DEG1 = c(DEG_genes_III_I$genes)
DEG2 = c(DEG_genes_TvC$genes)
DEG3 = c(DEG_genes_Paired$genes)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(DEG1, DEG2, DEG3), 
  category.names = c("Tumores", "Todos", "Pareados"),
  filename = "Venndiagram_DEG.png",
  output = TRUE,
  lwd = 2,  # Circles
  lty = 'blank',
  fill = myCol,
  cex = 1,  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("darkgreen", "darkred", "darkblue"),
  rotation = 1)



### Create a heatmap 
grouped_tumors <- arrange(Metadata, stage_group) %>% subset(grouped_tumors$Tissue == "Pancreatic tumor")
conds_tumors = grouped_tumors$stage_group
conds_tumors = as.factor(conds_tumors)
cond_colours_tumors = rainbow(length(unique(conds_tumors)), alpha = 1)[as.factor(conds_tumors)]
sig_genes_tumors <- subset(matrix_expression, subset = rownames(matrix_expression) %in% III_I$delabel)
sig_genes_tumors <- filter(sig_genes_tumors[, colnames(sig_genes_tumors) %in%
                                              grouped_tumors[grouped_tumors$Tissue == "Pancreatic tumor",]$geo_accession])
col_order <- c(grouped_tumors$geo_accession)
sig_genes_tumors <- sig_genes_tumors[ ,col_order]
                            
inter_tumors = unique(grouped_tumors$stage_group)
c = factor(conds_tumors, levels = inter_tumors) 

pdf("heatmap_Tumors_Group.pdf") 
heatmap(as.matrix(sig_genes_tumors), 
        col = hmcol, 
        cexRow = 0.8, 
        cexCol = 0.5,
        Colv = NA,
        ColSideColors = cond_colours_tumors)
legend("topleft", 
       levels(c), 
       fill = unique(cond_colours_tumors), 
       cex = 0.4)
heat_legend(as.matrix(sig_genes_tumors), 
            hmcol, 
            cex.lab = 0.4)
dev.off()

################################################################################
##                      4. FUNCTIONAL ENRICHMENT ANALYSIS                     ##
################################################################################

# a) KEGG

## Convert ID genes to Entrez
gene_id <- getBM(attributes = c("ensembl_gene_id",
                                "external_gene_name", 
                                "entrezgene_id"),
                 mart = useDataset("hsapiens_gene_ensembl", 
                                   useMart("ensembl")))

GeneID_TvC <- merge(TvC, gene_id[,c(2,3)], 
                    by.x = "genes", 
                    by.y = "external_gene_name")
GeneID_Paired <- merge(Paired_results, gene_id[,c(2,3)], 
                       by.x = "genes", 
                       by.y = "external_gene_name")
GeneID_II_I <- merge(II_I, gene_id[,c(2,3)], 
                     by.x = "genes", 
                     by.y = "external_gene_name")
GeneID_III_I <- merge(III_I, gene_id[,c(2,3)], 
                      by.x = "genes", 
                      by.y = "external_gene_name")
GeneID_IV_I <- merge(IV_I, gene_id[,c(2,3)], 
                     by.x = "genes", 
                     by.y = "external_gene_name")
GeneID_III_II <- merge(III_I, gene_id[,c(2,3)], 
                       by.x = "genes", 
                       by.y = "external_gene_name")
GeneID_IV_II <- merge(IV_II, gene_id[,c(2,3)], 
                      by.x = "genes", 
                      by.y = "external_gene_name")
GeneID_IV_III <- merge(IV_III, gene_id[,c(2,3)], 
                       by.x = "genes",
                       by.y = "external_gene_name")

#dim(GeneID_II_I) # 29670 obs. 10 variables
#dim(GeneID_III_I) # 29668 obs. 10 variables

## Remove rows with NA in the column entrezgene_id
GeneID_TvC <- GeneID_TvC %>% drop_na(entrezgene_id)
GeneID_Paired <- GeneID_Paired %>% drop_na(entrezgene_id)
GeneID_II_I <- GeneID_II_I %>% drop_na(entrezgene_id)
GeneID_III_I <- GeneID_III_I %>% drop_na(entrezgene_id)
GeneID_IV_I <- GeneID_IV_I %>% drop_na(entrezgene_id)
GeneID_III_II <- GeneID_III_II %>% drop_na(entrezgene_id)
GeneID_IV_II <- GeneID_IV_II %>% drop_na(entrezgene_id)
GeneID_IV_III <- GeneID_IV_III %>% drop_na(entrezgene_id)

#sum(is.na(GeneID_TvC$entrezgene_id)) # 0 (Example)
sum(is.na(GeneID_III_I$entrezgene_id))

## Generate Gene_lists for KEGG
gene_list_TvC = NULL  ##Tumor Vs Control

for(x in GeneID_TvC$entrezgene_id){
  gene_list_TvC = c(gene_list_TvC, GeneID_TvC[GeneID_TvC$entrezgene_id == x, "logFC"])
} 
names(gene_list_TvC) = GeneID_TvC$entrezgene_id

gene_list_TvC <- sort(gene_list_TvC, decreasing = TRUE)

head(gene_list_TvC)

gene_list_Paired = NULL  ##Paired samples

for(x in GeneID_Paired$entrezgene_id){
  gene_list_Paired = c(gene_list_Paired, GeneID_Paired[GeneID_Paired$entrezgene_id == x, "logFC"])
} 
names(gene_list_Paired) = GeneID_Paired$entrezgene_id

gene_list_Paired <- sort(gene_list_Paired, decreasing = TRUE)

head(gene_list_Paired)

gene_list_II_I = NULL  ##Tumor_II_I

for(x in GeneID_II_I$entrezgene_id){
  gene_list_II_I = c(gene_list_II_I, GeneID_II_I[GeneID_II_I$entrezgene_id == x, "logFC"])
} 

names(gene_list_II_I) = GeneID_II_I$entrezgene_id

gene_list_II_I <- sort(gene_list_II_I, decreasing = TRUE)

head(gene_list_II_I)


gene_list_III_I = NULL  ##Tumor_III_I

for(x in GeneID_III_I$entrezgene_id){
  gene_list_III_I = c(gene_list_III_I, GeneID_III_I[GeneID_III_I$entrezgene_id == x, "logFC"])
} 
names(gene_list_III_I) = GeneID_III_I$entrezgene_id

gene_list_III_I <- sort(gene_list_III_I, decreasing = TRUE)

head(gene_list_III_I)


gene_list_IV_I = NULL  ##Tumor_I_IV

for(x in GeneID_IV_I$entrezgene_id){
  gene_list_IV_I = c(gene_list_IV_I, GeneID_IV_I[GeneID_IV_I$entrezgene_id == x, "logFC"])
} 
names(gene_list_IV_I) = GeneID_IV_I$entrezgene_id

gene_list_IV_I <- sort(gene_list_IV_I, decreasing = TRUE)

head(gene_list_IV_I)


gene_list_III_II = NULL  ##Tumor_III_II

for(x in GeneID_III_II$entrezgene_id){
  gene_list_III_II = c(gene_list_III_II, GeneID_III_II[GeneID_III_II$entrezgene_id == x, "logFC"])
} 
names(gene_list_III_II) = GeneID_III_II$entrezgene_id

gene_list_III_II <- sort(gene_list_III_II, decreasing = TRUE)

head(gene_list_III_II)


gene_list_IV_II = NULL  ##Tumor_IV_II

for(x in GeneID_IV_II$entrezgene_id){
  gene_list_IV_II = c(gene_list_IV_II, GeneID_IV_II[GeneID_IV_II$entrezgene_id == x, "logFC"])
} 

names(gene_list_IV_II) = GeneID_IV_II$entrezgene_id

gene_list_IV_II <- sort(gene_list_IV_II, decreasing = TRUE)

head(gene_list_IV_II)


gene_list_IV_III = NULL  ##Tumor_IV_III

for(x in GeneID_IV_II$entrezgene_id){
  gene_list_IV_III = c(gene_list_IV_III, GeneID_IV_III[GeneID_IV_III$entrezgene_id == x, "logFC"])
} 
names(gene_list_IV_III) = GeneID_IV_III$entrezgene_id

gene_list_IV_III <- sort(gene_list_IV_III, decreasing = TRUE)

head(gene_list_IV_III)

## Calculate KEGG
TvC_KEGG <- gseKEGG(geneList = gene_list_TvC, 
                    organism = "hsa", 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")
head(TvC_KEGG)

Paired_KEGG <- gseKEGG(geneList = gene_list_Paired, 
                       organism = "hsa", 
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
head(Paired_KEGG)
dim(Paired_KEGG)

II_I_KEGG <- gseKEGG(geneList = gene_list_II_I, 
                     organism = "hsa", 
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")
head(II_I_KEGG)

III_I_KEGG <- gseKEGG(geneList = gene_list_III_I, 
                      organism = "hsa", 
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")
head(III_I_KEGG)

IV_I_KEGG <- gseKEGG(geneList = gene_list_IV_I, 
                     organism = "hsa", 
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")
head(IV_I_KEGG)

III_II_KEGG <- gseKEGG(geneList = gene_list_III_II, 
                       organism = "hsa", 
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
head(III_II_KEGG)

IV_II_KEGG <- gseKEGG(geneList = gene_list_IV_II, 
                      organism = "hsa", 
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")
head(IV_II_KEGG)

IV_III_KEGG <- gseKEGG(geneList = gene_list_IV_III, 
                       organism = "hsa", 
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
head(IV_III_KEGG)

head(gene_list_III_I)

## dotplots KEGG
dotplot(TvC_KEGG, #ALL
        showCategory = 10,
        split =".sign",
        title = "KEGG Todos",
        font.size = 7,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

dotplot(Paired_KEGG, #Paired
        showCategory = 10,
        split =".sign",
        title = "KEGG Pareados",
        font.size = 7,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

dotplot(III_I_KEGG, #III_I
        showCategory = 10,
        split =".sign",
        title = "KEGG Tumores",
        font.size = 10,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "Tumors_KEGG.pdf")
dev.off()


# b) GO (Gene Onthology)

## Generate Gene_lists for GO with symbols
gene_list_TvC_GO = NULL  ##Tumor Vs Control

for(x in GeneID_TvC$genes){
  gene_list_TvC_GO = c(gene_list_TvC_GO, GeneID_TvC[GeneID_TvC$genes == x, "logFC"])
} 
names(gene_list_TvC_GO) = GeneID_TvC$genes

gene_list_TvC_GO <- sort(gene_list_TvC_GO, decreasing = TRUE)


gene_list_Paired_GO = NULL  ##Paired samples

for(x in GeneID_Paired$genes){
  gene_list_Paired_GO = c(gene_list_Paired_GO, GeneID_Paired[GeneID_Paired$genes == x, "logFC"])
} 
names(gene_list_Paired_GO) = GeneID_Paired$genes

gene_list_Paired_GO <- sort(gene_list_Paired_GO, decreasing = TRUE)


gene_list_II_I_GO = NULL  ##Tumor_II_I

for(x in GeneID_II_I$genes){
  gene_list_II_I_GO = c(gene_list_II_I_GO, GeneID_II_I[GeneID_II_I$genes == x, "logFC"])
} 
names(gene_list_II_I_GO) = GeneID_II_I$genes

gene_list_II_I_GO <- sort(gene_list_II_I_GO, decreasing = TRUE)

head(gene_list_II_I_GO)


gene_list_III_I_GO = NULL  ##Tumor_III_I

for(x in GeneID_III_I$genes){
  gene_list_III_I_GO = c(gene_list_III_I_GO, GeneID_III_I[GeneID_III_I$genes == x, "logFC"])
} 
names(gene_list_III_I_GO) = GeneID_III_I$genes

gene_list_III_I_GO <- sort(gene_list_III_I_GO, decreasing = TRUE)

head(gene_list_III_I_GO)


gene_list_IV_I_GO = NULL  ##Tumor_IV_I

for(x in GeneID_IV_I$genes){
  gene_list_IV_I_GO = c(gene_list_IV_I_GO, GeneID_IV_I[GeneID_IV_I$genes == x, "logFC"])
} 
names(gene_list_IV_I_GO) = GeneID_IV_I$genes

gene_list_IV_I_GO <- sort(gene_list_IV_I_GO, decreasing = TRUE)

head(gene_list_IV_I_GO)


gene_list_III_II_GO = NULL  ##Tumor_III_II

for(x in GeneID_III_II$genes){
  gene_list_III_II_GO = c(gene_list_III_II_GO, GeneID_III_II[GeneID_III_II$genes == x, "logFC"])
} 
names(gene_list_III_II_GO) = GeneID_III_II$genes

gene_list_III_II_GO <- sort(gene_list_III_II_GO, decreasing = TRUE)

head(gene_list_III_II_GO)


gene_list_IV_II_GO = NULL  ##Tumor_IV_II

for(x in GeneID_IV_II$genes){
  gene_list_IV_II_GO = c(gene_list_IV_II_GO, GeneID_IV_II[GeneID_IV_II$genes == x, "logFC"])
} 
names(gene_list_IV_II_GO) = GeneID_IV_II$genes

gene_list_IV_II_GO <- sort(gene_list_IV_II_GO, decreasing = TRUE)

head(gene_list_IV_II_GO)


gene_list_IV_III_GO = NULL  ##Tumor_IV_III

for(x in GeneID_IV_III$genes){
  gene_list_IV_III_GO = c(gene_list_IV_III_GO, GeneID_IV_III[GeneID_IV_III$genes == x, "logFC"])
} 
names(gene_list_IV_III_GO) = GeneID_IV_III$genes

gene_list_IV_III_GO <- sort(gene_list_IV_III_GO, decreasing = TRUE)

head(gene_list_IV_III_GO)

## Calculate ontologies GO
TvC_GO_MF <- gseGO(geneList = gene_list_TvC_GO, ##Tumor vs Control (TvC)
                   ont = "MF",
                   OrgDb = "org.Hs.eg.db",    
                   keyType = "SYMBOL",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 120)
#head(TvC_GO_MF)

TvC_GO_BP <- gseGO(geneList = gene_list_TvC_GO,
                   ont = "BP",
                   OrgDb = "org.Hs.eg.db",    
                   keyType = "SYMBOL",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 120)
#head(TvC_GO_BP)

TvC_GO_CC <- gseGO(geneList = gene_list_TvC_GO,
                   ont = "CC",
                   OrgDb = "org.Hs.eg.db",    
                   keyType = "SYMBOL",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 120)
#head(TvC_GO_CC)


Paired_GO_MF <- gseGO(geneList = gene_list_Paired_GO,  ## Paired samples
                      ont = "MF",
                      OrgDb = "org.Hs.eg.db",    
                      keyType = "SYMBOL",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      maxGSSize = 120)
#head(Paired_GO_MF)

Paired_GO_BP <- gseGO(geneList = gene_list_Paired_GO,
                      ont = "BP",
                      OrgDb = "org.Hs.eg.db",    
                      keyType = "SYMBOL",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      maxGSSize = 120)
#head(Paired_GO_BP)

Paired_GO_CC <- gseGO(geneList = gene_list_Paired_GO,
                      ont = "CC",
                      OrgDb = "org.Hs.eg.db",    
                      keyType = "SYMBOL",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      maxGSSize = 120)
#head(Paired_GO_CC)


II_I_GO_MF <- gseGO(geneList = gene_list_II_I_GO, #Tumor II_I
                    ont = "MF",
                    OrgDb = "org.Hs.eg.db",    
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 120)
#head(II_I_GO_MF)


II_I_GO_BP <- gseGO(geneList = gene_list_II_I_GO,
                    ont = "BP",
                    OrgDb = "org.Hs.eg.db",    
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 120)
#head(II_I_GO_BP)

II_I_GO_CC <- gseGO(geneList = gene_list_II_I_GO,
                    ont = "CC",
                    OrgDb = "org.Hs.eg.db",    
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 120)
#head(II_I_GO_CC)


III_I_GO_MF <- gseGO(geneList = gene_list_III_I_GO, #Tumor III_I
                    ont = "MF",
                    OrgDb = "org.Hs.eg.db",    
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 120)
#head(III_I_GO_MF)

III_I_GO_BP <- gseGO(geneList = gene_list_III_I_GO,
                    ont = "BP",
                    OrgDb = "org.Hs.eg.db",    
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 120)
#head(III_I_GO_BP)

III_I_GO_CC <- gseGO(geneList = gene_list_III_I_GO,
                    ont = "CC",
                    OrgDb = "org.Hs.eg.db",    
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 120)
#head(III_I_GO_CC)


## Dotplot GO
dotplot(TvC_GO_MF, ##ALL: Tumours vs Control (TvC)
        showCategory = 10,
        split =".sign",
        title = "Enriched GO Molecular Function",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 6))
pdf(file = "GO_TvC_MF.pdf")
dev.off()

dotplot(TvC_GO_BP,
        showCategory = 10,
        split =".sign",
        title = "Enriched GO Biological process",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "GO_TvC_BP.pdf")
dev.off()

dotplot(TvC_GO_CC,
        showCategory = 10,
        split =".sign",
        title = "Enriched GO cellular component",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "GO_TvC_CC.pdf")
dev.off()


dotplot(Paired_GO_MF, ##Paired
        showCategory = 10,
        split =".sign",
        title = "Enriched GO Molecular Function",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 6))
pdf(file = "GOpairedMF.pdf")
dev.off()

dotplot(Paired_GO_BP,
        showCategory = 10,
        split =".sign",
        title = "Enriched GO Biological process",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "GOpairedBP.pdf")
dev.off()

dotplot(Paired_GO_CC,
        showCategory = 10,
        split =".sign",
        title = "Enriched GO cellular component",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "GOpairedCC.pdf")
dev.off()

dotplot(III_I_GO_MF, ##Tumors
        showCategory = 10,
        split =".sign",
        title = "Enriched GO Molecular Function",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 6))
pdf(file = "GO-III_I-MF.pdf")
dev.off()

dotplot(III_I_GO_BP,
        showCategory = 10,
        split =".sign",
        title = "Enriched GO Biological process",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "GO-III_I-BP.pdf")
dev.off()

dotplot(III_I_GO_CC,
        showCategory = 10,
        split =".sign",
        title = "Enriched GO cellular component",
        font.size = 8,
        orderBy = "x") +
  facet_grid(.~.sign) +
  scale_color_gradient(low = "#E74C3C", high ="#438FE6") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
pdf(file = "GO-III_I-CC.pdf")
dev.off()








