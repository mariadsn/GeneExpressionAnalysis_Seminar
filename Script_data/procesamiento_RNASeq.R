## DATOS DE RNA-Seq

# instalacion de Bioconductor
install.packages("BiocManager") # gestiona paquetes de bioconductor
BiocManager::install("maftools")
BiocManager::install("TCGAbiolinks", forced=TRUE)

library(BiocManager)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

library(pheatmap)
library(maftools)

library(limma)

library(DT)
library(NOISeq)
library(Cairo)
library(scales)
library(ComplexHeatmap)
library(EDASeq)


## RECOPILACION
# descarga de datos TCGA-SKCM
# obtener lista de los proyectos disponibles en GDC
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-SKCM') # resumen del proyecto

# Lanza la consulta para recuperar todos los datos que satisfagan los criterios de búsqueda (breast cancer data)
query_TCGA <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification") # aplicacion de filtros
output_query_TCGA <- getResults(query_TCGA)
# lanzar consulta para: gene expression data; estrategia experimental: RNA-Seq, flujo de trabajo de analisis: STAR-counts, acceso: open
query_TCGA <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       access = "open",
                       barcode = c("TCGA-BH-A18M-11A-33R-A12D-07", "TCGA-BH-A18R-11A-42R-A12D-07", 
                                   "TCGA-BH-A1FE-11B-14R-A13Q-07", "TCGA-BH-A0BZ-11A-61R-A12P-07",
                                   "TCGA-BH-A0DH-11A-31R-A089-07", "TCGA-BH-A1F0-11B-23R-A137-07",
                                   "TCGA-PL-A8LV-01A-21R-A41B-07", "TCGA-BH-A0BC-01A-22R-A084-07",
                                   "TCGA-AR-A1AX-01A-11R-A12P-07", "TCGA-AC-A2FO-01A-11R-A180-07",
                                   "TCGA-AQ-A0Y5-01A-11R-A14M-07", "TCGA-AC-A3EH-01A-22R-A22K-07"))
getResults(query_TCGA)
# descarga de los datos
GDCdownload(query_TCGA)

#preprocesamiento de los datos descargados (counts)
SKCM.counts <- GDCprepare(query = query_TCGA,
                          summarizedExperiment = TRUE)
rm(query_TCGA)

# Matriz de expresión
counts_data <- assay(SKCM.counts)

colData <- read.delim("id_infoStudy.csv", stringsAsFactors = TRUE, sep = ',')

# convert counts to DGEList object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ sample_type)
# prefiltering: removing rows with low gene counts
# keeping rown that have at least 10 reads total
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# set the factor level
dds$sample_type <- relevel(dds$sample_type, ref = 'Solid Tissue Normal')

# Run DESeq
dds <- DESeq(dds) # datos para el analisis diferencial
res <- results(dds)
summary(res)
countdata <- counts(dds) # para obtener el log2
countdata_log2 <- log2(countdata +1)

# explore results
fold.change.Breast <- res$log2FoldChange # para obtener FC
adj.Pval.Breast <- res$padj # pvalor ajustado
genes.ids.Breast <- rownames(res)
# metodo combinado para los genes activados y reprimidos
activated.genes.breast.1 <- genes.ids.Breast[fold.change.Breast > 2 & adj.Pval.Breast < 0.005]
repressed.genes.breast.1 <- genes.ids.Breast[fold.change.Breast < -2 & adj.Pval.Breast < 0.005]
length(activated.genes.breast.1)
length(repressed.genes.breast.1)
# visualizacion de los datos volcano plot
names(fold.change.Breast) <- genes.ids.Breast
log.padj.breast <- -log10(adj.Pval.Breast)
names(log.padj.breast) <- genes.ids.Breast
plot(fold.change.Breast,log.padj.breast, ylab="-log10(p value)", xlab= "log2 fold change", pch=19, cex=0.5,
     col = "grey", xlim=c(-6,6))
points(fold.change.Breast[activated.genes.breast.1],log.padj.breast[activated.genes.breast.1], 
       pch = 19, cex = 0.5, col = "red")
points(fold.change.Breast[repressed.genes.breast.1],log.padj.breast[repressed.genes.breast.1], 
       pch = 19, cex = 0.5, col = "blue")

# obtencion de los datos log2
breast.all.DEG <- union(activated.genes.breast.1, repressed.genes.breast.1)
# cambio nombres columnas
colnames(countdata_log2) <-sampleinfo$name_id
normal.breast.DEG.table <- countdata_log2[, c("Normal_breast_1", "Normal_breast_2", 
                                                "Normal_breast_3", "Normal_breast_4", 
                                                "Normal_breast_5", "Normal_breast_6")]
normal.breast.DEG.table <- normal.breast.DEG.table[rownames(normal.breast.DEG.table) %in% breast.all.DEG,]
normal.breast.DEG.table <- cbind(attr_name = rownames(normal.breast.DEG.table), normal.breast.DEG.table)

tumor.breast.DEG.table <- countdata_log2[, c("Tumor_breast_1", "Tumor_breast_2", 
                                               "Tumor_breast_3", "Tumor_breast_4", 
                                               "Tumor_breast_5", "Tumor_breast_6")]
tumor.breast.DEG.table <- tumor.breast.DEG.table[rownames(tumor.breast.DEG.table) %in% breast.all.DEG,]
tumor.breast.DEG.table <- cbind(attr_name = rownames(tumor.breast.DEG.table), tumor.breast.DEG.table)

write.table(normal.breast.DEG.table, file="../NormalBreast_DEGs_table.csv", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tumor.breast.DEG.table, file="../TumorBreast_DEGs_table.csv", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)




