## Procesamiento datos Microarray

# instalacion de Bioconductor
install.packages("BiocManager") # gestiona paquetes de bioconductor
# instalacion de paquetes a partir de BiocManager
BiocManager::install("affy", force = TRUE)
BiocManager::install("affyPLM", force = TRUE)
BiocManager::install("arrayQualityMetrics", force = TRUE)
BiocManager::install("limma", force = TRUE)
BiocManager::install("annaffy", force = TRUE)
BiocManager::install("hgu133plus2.db", force = TRUE)

# fijar el directorio de trabajo a la carpeta con el script de R
# la carpeta que contiene los dataset tambien debe de incluir este script

# iniciacion de los paquetes
library("affy")
library("affyPLM")
library("arrayQualityMetrics")
library("limma")
library("annaffy")
library("hgu133plus2.db")


# carga de ficheros con los datos en bruto (formato CEL)
SDRF <- read.delim("../data_info.csv", sep = ",")
rownames(SDRF) <- SDRF$Ids
SDRF <- AnnotatedDataFrame(SDRF)

microarray.raw.data <- ReadAffy(filenames = SDRF$Array.Data.File, verbose=TRUE, phenoData = SDRF) 
microarray.raw.data

### QUALITY CONTROL
# imagenes de escaneres (daños fisicos)
image(microarray.raw.data[,1], col=rainbow(100))
image(microarray.raw.data[,2], col=rainbow(100))
image(microarray.raw.data[,3], col=rainbow(100))
image(microarray.raw.data[,4], col=rainbow(100))
image(microarray.raw.data[,5], col=rainbow(100))
image(microarray.raw.data[,6], col=rainbow(100))
image(microarray.raw.data[,7], col=rainbow(100))
image(microarray.raw.data[,8], col=rainbow(100))
image(microarray.raw.data[,9], col=rainbow(100))
image(microarray.raw.data[,10], col=rainbow(100))
image(microarray.raw.data[,11], col=rainbow(100))
image(microarray.raw.data[,12], col=rainbow(100))

# metricas de calidad
arrayQualityMetrics(expressionset = microarray.raw.data,
                    outdir = "../dataArrayQualityMetrics_whAffy",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = "samples")

### DATA PREPROCESSING
# reducir el ruido por variabilidad experimental pero manteniendo la variabilidad biologica 
# antes de eliminar el ruido
boxplot(microarray.raw.data, col = rainbow(12), las=2, ylab="Fluorescence")
# Robust multiarray average (rma) realiza la correccion de la fluorescencia de fondo, 
# normaliza y calcula los niveles de expresion con transformacion de los datos en log2
microarray.processed.data <-rma(microarray.raw.data)
# despues de eliminar el ruido y normalizar (datos comparables entre si)
boxplot(microarray.processed.data, col = rainbow(12), las=2, ylab="Fluorescence A.U.")
# extraccion de la matriz de interes
expression.level <- exprs(microarray.processed.data)
head(expression.level)
dim(expression.level)
# modificar nombre de las columnas
sampleID <- c("normal_breast_1","tumor_breast_1","normal_breast_2","tumor_breast_2",
              "normal_breast_3","tumor_breast_3","normal_breast_4","tumor_breast_4",
              "normal_breast_5","tumor_breast_5","normal_breast_6","tumor_breast_6")
colnames(expression.level) <- sampleID
head(expression.level)
# calculo de la media entre condiciones (facilita el analisis comparativo)
normal.breast <- (expression.level[,"normal_breast_1"] + expression.level[,"normal_breast_2"] + 
                    expression.level[,"normal_breast_3"] + expression.level[,"normal_breast_4"] +
                    expression.level[,"normal_breast_5"] + expression.level[,"normal_breast_6"])/6
tumor.breast <- (expression.level[,"tumor_breast_1"] + expression.level[,"tumor_breast_2"] + 
                   expression.level[,"tumor_breast_3"] + expression.level[,"tumor_breast_4"] +
                   expression.level[,"tumor_breast_5"] + expression.level[,"tumor_breast_6"])/6
# generacion de matriz que guarda los datos de la media
mean.expression <- matrix(c(normal.breast, tumor.breast), ncol=2)
conditions.id <- c("normal_breast", "tumor_breast")
rownames(mean.expression) <- names(normal.breast)
colnames(mean.expression) <- conditions.id
head(mean.expression)
# grafico de dispersion para un analisis exploratorio preliminar
plot(tumor.breast,normal.breast, xlab="Normal stromal breast", ylab= "Tumor stromal breast", pch=19, cex=0.5)

### INFERENTIAL STATISTICS
# generacion de la matriz que contiene el diseño experimental
experimental.desing <- model.matrix(~ -1+factor(c(1,2,1,2,1,2,1,2,1,2,1,2)))
colnames(experimental.desing) <- c("normal_breast", "tumor_breast")
# Seleccion de genes que se expresan de forma diferencial (limma):
linear.fit <- lmFit(expression.level, experimental.desing) # calcula la media de cada condicion
contrast.matrix <-makeContrasts(tumor_breast- normal_breast, 
                                levels = c("normal_breast","tumor_breast"))  # especificacion de contrastes
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix) # calculo del Fold-Change
contrast.results <- eBayes(contrast.linear.fit) # calculo p-valor/q-valor
# extraer DEGs -- breast
Breast.results <- topTable(contrast.results, number = 54675, coef = 1, sort.by = "logFC")
head(Breast.results)
fold.change.Breast <- Breast.results$logFC # para obtener FC
adj.Pval.Breast <- Breast.results$adj.P.Val # pvalor ajustado
genes.ids.Breast <- rownames(Breast.results)
# metodo combinado para los genes activados y reprimidos
activated.genes.breast.1 <- genes.ids.Breast[fold.change.Breast > 1 & adj.Pval.Breast < 0.05]
repressed.genes.breast.1 <- genes.ids.Breast[fold.change.Breast < -1 & adj.Pval.Breast < 0.05]
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

### EXPLORATORY ANALYSIS
# obtencion del nombre de los genes
activated.genes.breast.table <- aafTableAnn(activated.genes.breast.1, "hgu133plus2.db", aaf.handler())
saveText(activated.genes.breast.table, file="../activatedGenes_breast_table.txt")

repressed.genes.breast.table <- aafTableAnn(repressed.genes.breast.1, "hgu133plus2.db", aaf.handler())
saveText(repressed.genes.breast.table, file="../repressedGenes_breast_table.txt")

# obtencion de los datos log2
breast.all.DEG <- union(activated.genes.breast.1, repressed.genes.breast.1)
normal.breast.DEG.table <- expression.level[, c("normal_breast_1", "normal_breast_2", 
                                                "normal_breast_3", "normal_breast_4", 
                                                "normal_breast_5", "normal_breast_6")]
normal.breast.DEG.table <- normal.breast.DEG.table[rownames(normal.breast.DEG.table) %in% breast.all.DEG,]
normal.breast.DEG.table <- cbind(attr_name = rownames(normal.breast.DEG.table), normal.breast.DEG.table)

tumor.breast.DEG.table <- expression.level[, c("tumor_breast_1", "tumor_breast_2", 
                                               "tumor_breast_3", "tumor_breast_4", 
                                               "tumor_breast_5", "tumor_breast_6")]
tumor.breast.DEG.table <- tumor.breast.DEG.table[rownames(tumor.breast.DEG.table) %in% breast.all.DEG,]
tumor.breast.DEG.table <- cbind(attr_name = rownames(tumor.breast.DEG.table), tumor.breast.DEG.table)

write.table(normal.breast.DEG.table, file="../normalBreast_DEGs_table.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tumor.breast.DEG.table, file="../tumorBreast_DEGs_table.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
