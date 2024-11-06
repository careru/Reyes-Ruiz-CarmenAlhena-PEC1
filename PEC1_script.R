rm(list = ls())

getwd()
setwd("C:/Users/USUARIO/Dropbox/MASTER_UOC/4. SEMESTRE 2025-2025/Análisis de datos ómicos (repitiendo)/PEC1")
getwd()


#Cargar el los datos en txt
library(readr)
data <- read.table("ST000002_AN000002_clean.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Extraer y preparar los datos
sample_names <- as.character(data[1, -1])  # Crea un vector con nombres de las muestras de la fila 1
groups <- as.character(data[2, -1])        # Crea un vector con nombres de los grupos de la fila 2
row_names <- data[-c(1,2), 1]              # Almacenar los nombres de los metabolitos, primera columna
counts <- data[-c(1,2), -1]                # Elimina las dos primeras filas y la primera columna para quedarse solo con los datos
counts <- apply(counts, 2, as.numeric)     # Convertir los datos a formato numérico
rownames(counts) <- row_names              # Asigna los nombres de las filas a la matriz de conteos
colnames(counts) <- sample_names           # Asigna los nombres de las columnas a la matriz de conteos

# Crear el DataFrame de metadatos para las muestras
col_data <- data.frame(Group = groups, row.names = sample_names)

# Crear el objeto SummarizedExperiment (se) con el DataFrame y los datos numéricos
library(SummarizedExperiment)
se <- SummarizedExperiment(
  assays = list(counts = counts),
  colData = col_data
)


#Exploarar la clase de datos
se 
colData(se)  
se$Group  

rownames(se) 

assay(se)



# Análisis de PCA

# PCA y dataframe de resultados
pca <- prcomp(t(assay(se)), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Group <- colData(se)$Group

# Gráfico de PCA
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA metabolomics of intestinal tissue",
       x = "PC1", y = "PC2") +
  scale_color_manual(values = c("Before" = "blue", "After" = "orange"))


# Análisis diferencial

# Instalar y cargar limma
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("limma")
}
library(limma)
# Crear la matriz (m1) y ajuste (fit)
m1 <- model.matrix(~ colData(se)$Group)
fit <- lmFit(assay(se), m1 )
fit <- eBayes(fit)
#Resultados significativos
results <- topTable(fit, coef = 2)
print(results)

