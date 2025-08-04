

#  Setup personal library path  
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



#   Install TCGAbiolinks and other required Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), ask = FALSE)

# Load the libraries
library(TCGAbiolinks)
library(SummarizedExperiment)


#   Query and download TCGA BRCA gene expression data 

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query)

# Load the data into R as a Summarized Experiment object 
data <- GDCprepare(query)



#  Cell 4: Explore sample metadata and count matrix

colData(data)

# View count matrix (gene expression)
assay(data)[1:5, 1:5]  # Shows first 5 genes across 5 samples


#   Install and load DESeq2, create DESeqDataSet and filter 


if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)


dds <- DESeqDataSet(data, design = ~ sample_type)

colnames(colData(data))


dds <- dds[rowSums(counts(dds)) > 10, ]

#   Run differential expression analysis with DESeq2 

dds <- DESeq(dds)

# Extract results
res <- results(dds)

# View summary
summary(res)

# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# View top results
head(resOrdered)



#  Convert to data frame and write to CSV
resDF <- as.data.frame(resOrdered)
write.csv(resDF, file = "DESeq2_results.csv")

#  Plot MA plot and volcano plot 

plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))

# Basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot", col=ifelse(padj<0.05, "red", "black")))



# Install and load annotation and enrichment packages 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "enrichplot", "pheatmap", "DESeq2", "TCGAbiolinks"), ask = FALSE)

library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(DESeq2)
library(TCGAbiolinks)



# Assuming we have  already have resDF from DESeq2 results ---

#  Annotate DESeq2 results with gene symbols and names
clean_ids <- gsub("\\.\\d+$", "", rownames(resDF))

# Add gene symbols and names to DESeq2 results
resDF$symbol <- mapIds(
  org.Hs.eg.db,
  keys = clean_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

resDF$genename <- mapIds(
  org.Hs.eg.db,
  keys = clean_ids,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(resDF[, c("symbol", "genename")])



# Prepare significant genes for enrichment analysis 

sigGenes <- rownames(resDF)[which(!is.na(resDF$padj) & resDF$padj < 0.05)]

# Clean Ensembl IDs for enrichment analysis
sigGenes_clean <- gsub("\\.\\d+$", "", sigGenes)

# Convert Ensembl IDs to Entrez IDs for clusterProfiler
entrezIDs <- bitr(
  sigGenes_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#  GO enrichment analysis (Biological Process) ---
ego <- enrichGO(
  gene = entrezIDs$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE  # Show gene symbols in output
)

# Visualize top 20 enriched GO terms
dotplot(ego, showCategory = 20) + ggtitle("Top 20 Enriched GO Biological Processes")

# KEGG Pathway enrichment analysis
ekegg <- enrichKEGG(
  gene = entrezIDs$ENTREZID,
  organism = "hsa",  # Human
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Visualize top 10 KEGG pathways
barplot(ekegg, showCategory = 10) + ggtitle("Top 10 Enriched KEGG Pathways")

# Heatmap of top 50 differentially expressed genes
norm_counts <- counts(dds, normalized = TRUE)
topGenes <- head(rownames(resOrdered), 50)
mat <- norm_counts[topGenes, ]

# Scale rows (genes) for visualization
mat_scaled <- t(scale(t(mat)))

pheatmap(
  mat_scaled,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Genes Heatmap"
)

# PCA plot to visualize tumor vs normal separation
vsd <- vst(dds, blind = FALSE)  # variance stabilizing transformation

plotPCA(vsd, intgroup = "sample_type") + ggtitle("PCA: Tumor vs Normal Samples")



# Save annotated DE results
write.csv(resDF, "DESeq2_Annotated_Results.csv")

# Save GO and KEGG enrichment results
write.csv(ego@result, "GO_Enrichment_Results.csv")
write.csv(ekegg@result, "KEGG_Enrichment_Results.csv")





# 1. Create variance stabilized data (vsd) from dds
vsd <- vst(dds, blind = FALSE)

# 2. Extract expression matrix and labels (use correct label column name)
expr_matrix <- assay(vsd)
labels <- colData(vsd)$sample_type  # this is the group info: tumor vs normal

# 3. Prepare dataframe for ML (samples as rows)
expr_df <- as.data.frame(t(expr_matrix))
expr_df$condition <- as.factor(labels)  # add labels as factor column

# 4. Load caret and split dataset into train/test sets
library(caret)
set.seed(123)
trainIndex <- createDataPartition(expr_df$condition, p = 0.8, list = FALSE)

trainData <- expr_df[trainIndex, ]
testData <- expr_df[-trainIndex, ]



# Calculate variance of each gene across samples
gene_vars <- apply(trainData[, -ncol(trainData)], 2, var)

# Select top 500 most variable genes
topGenes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# Subset train and test data to only these genes + condition
trainData_sub <- trainData[, c(topGenes, "condition")]
testData_sub <- testData[, c(topGenes, "condition")]

install.packages("randomForest")
library(randomForest)




# Train random forest on reduced features
set.seed(123)
rf_model <- randomForest(
  condition ~ .,
  data = trainData_sub,
  importance = TRUE,
  ntree = 500
)

# Predict and evaluate
predictions <- predict(rf_model, newdata = testData_sub)
confusionMatrix(predictions, testData_sub$condition)

# View model summary
print(rf_model)

# Predict on test data
predictions <- predict(rf_model, newdata = testData)

# Confusion matrix to evaluate performance
confusionMatrix(predictions, testData$condition)




# Plot top 20 important features
varImpPlot(rf_model, n.var = 20, main = "Top 20 Important Genes")



library(pROC)

# Get probabilities
probs <- predict(rf_model, newdata = testData_sub, type = "prob")[,2]

# Plot ROC
roc_obj <- roc(testData_sub$condition, probs)
plot(roc_obj, main = "ROC Curve")
auc(roc_obj)



saveRDS(rf_model, file = "random_forest_model.rds")
# Later you can load it using:
# model <- readRDS("random_forest_model.rds")








# Export expression values + labels
write.csv(expr_df, "expression_data.csv", row.names = TRUE)

# Export random forest predictions
results <- data.frame(
  Sample = rownames(testData),
  Actual = testData$condition,
  Predicted = predictions
)
write.csv(results, "rf_predictions.csv", row.names = FALSE)

# Export feature importance
imp <- importance(rf_model)
imp_df <- data.frame(Gene = rownames(imp), MeanDecreaseGini = imp[, "MeanDecreaseGini"])
write.csv(imp_df, "gene_importance.csv", row.names = FALSE)

# Export GO / KEGG enrichment results
write.csv(ego@result, "GO_Enrichment.csv", row.names = FALSE)




# Save annotated DESeq2 results
write.csv(resDF, file = "DESeq2_Annotated_Results.csv", row.names = TRUE)


getwd()


write.csv(resDF, file = "C:/Users/Ganga Supriya/OneDrive/Documents/DESeq2_Annotated_Results.csv", row.names = TRUE)

save.image(file = "my_workspace.RData")



# Set personal library path if needed (especially if permission issues occur)
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install TCGAbiolinks and other required Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), ask = FALSE)

# Load the libraries
library(TCGAbiolinks)
library(SummarizedExperiment)



query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query)

# Load the data into R as a Summarized Experiment object 
data <- GDCprepare(query)



# View sample metadata
colData(data)

# View count matrix (gene expression)
assay(data)[1:5, 1:5]  # Shows first 5 genes across 5 samples



if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)


dds <- DESeqDataSet(data, design = ~ sample_type)

colnames(colData(data))


dds <- dds[rowSums(counts(dds)) > 10, ]


dds <- DESeq(dds)


# Extract results
res <- results(dds)

# View summary
summary(res)



# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# View top results
head(resOrdered)



# Convert to data frame and write to CSV
resDF <- as.data.frame(resOrdered)
write.csv(resDF, file = "DESeq2_results.csv")



plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))


# Basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano Plot", col=ifelse(padj<0.05, "red", "black")))



# Load necessary libraries (install once if needed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "enrichplot", "pheatmap", "DESeq2", "TCGAbiolinks"), ask = FALSE)

library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)
library(DESeq2)
library(TCGAbiolinks)

# --- Assume you already have resDF from DESeq2 results ---

# Remove version numbers from Ensembl IDs in results
clean_ids <- gsub("\\.\\d+$", "", rownames(resDF))

# Add gene symbols and names to DESeq2 results
resDF$symbol <- mapIds(
  org.Hs.eg.db,
  keys = clean_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

resDF$genename <- mapIds(
  org.Hs.eg.db,
  keys = clean_ids,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(resDF[, c("symbol", "genename")])

# Filter significant genes with adjusted p-value < 0.05
sigGenes <- rownames(resDF)[which(!is.na(resDF$padj) & resDF$padj < 0.05)]

# Clean Ensembl IDs for enrichment analysis
sigGenes_clean <- gsub("\\.\\d+$", "", sigGenes)

# Convert Ensembl IDs to Entrez IDs for clusterProfiler
entrezIDs <- bitr(
  sigGenes_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Gene Ontology (GO) enrichment analysis - Biological Process
ego <- enrichGO(
  gene = entrezIDs$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE  # Show gene symbols in output
)

# Visualize top 20 enriched GO terms
dotplot(ego, showCategory = 20) + ggtitle("Top 20 Enriched GO Biological Processes")

# KEGG Pathway enrichment analysis
ekegg <- enrichKEGG(
  gene = entrezIDs$ENTREZID,
  organism = "hsa",  # Human
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Visualize top 10 KEGG pathways
barplot(ekegg, showCategory = 10) + ggtitle("Top 10 Enriched KEGG Pathways")

# Heatmap of top 50 differentially expressed genes
norm_counts <- counts(dds, normalized = TRUE)
topGenes <- head(rownames(resOrdered), 50)
mat <- norm_counts[topGenes, ]

# Scale rows (genes) for visualization
mat_scaled <- t(scale(t(mat)))

pheatmap(
  mat_scaled,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 Differentially Expressed Genes Heatmap"
)

# PCA plot to visualize tumor vs normal separation
vsd <- vst(dds, blind = FALSE)  # variance stabilizing transformation

plotPCA(vsd, intgroup = "sample_type") + ggtitle("PCA: Tumor vs Normal Samples")




# Prepare data for machine learning analysis


#  Create variance stabilized data (vsd) from dds
vsd <- vst(dds, blind = FALSE)

#  Extract expression matrix and labels (use correct label column name)
expr_matrix <- assay(vsd)
labels <- colData(vsd)$sample_type  # this is the group info: tumor vs normal

#  Prepare dataframe for ML (samples as rows)
expr_df <- as.data.frame(t(expr_matrix))
expr_df$condition <- as.factor(labels)  # add labels as factor column

#  Load caret and split dataset into train/test sets
library(caret)
set.seed(123)
trainIndex <- createDataPartition(expr_df$condition, p = 0.8, list = FALSE)

trainData <- expr_df[trainIndex, ]
testData <- expr_df[-trainIndex, ]


# Calculate variance of each gene across samples
gene_vars <- apply(trainData[, -ncol(trainData)], 2, var)

# Select top 500 most variable genes
topGenes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# Subset train and test data to only these genes + condition
trainData_sub <- trainData[, c(topGenes, "condition")]
testData_sub <- testData[, c(topGenes, "condition")]

install.packages("randomForest")
library(randomForest)

# Train random forest on reduced features
set.seed(123)
rf_model <- randomForest(
  condition ~ .,
  data = trainData_sub,
  importance = TRUE,
  ntree = 500
)

# Predict and evaluate
predictions <- predict(rf_model, newdata = testData_sub)
confusionMatrix(predictions, testData_sub$condition)

# View model summary
print(rf_model)

# Predict on test data
predictions <- predict(rf_model, newdata = testData)

# Confusion matrix to evaluate performance
confusionMatrix(predictions, testData$condition)


# Plot top 20 important features
varImpPlot(rf_model, n.var = 20, main = "Top 20 Important Genes")

library(pROC)

# Get probabilities
probs <- predict(rf_model, newdata = testData_sub, type = "prob")[,2]

# Plot ROC
roc_obj <- roc(testData_sub$condition, probs)
plot(roc_obj, main = "ROC Curve")
auc(roc_obj)

saveRDS(rf_model, file = "random_forest_model.rds")


# Export expression values + labels
write.csv(expr_df, "expression_data.csv", row.names = TRUE)

# Export random forest predictions
results <- data.frame(
  Sample = rownames(testData),
  Actual = testData$condition,
  Predicted = predictions
)
write.csv(results, "rf_predictions.csv", row.names = FALSE)

# Export feature importance
imp <- importance(rf_model)
imp_df <- data.frame(Gene = rownames(imp), MeanDecreaseGini = imp[, "MeanDecreaseGini"])
write.csv(imp_df, "gene_importance.csv", row.names = FALSE)

# Export GO / KEGG enrichment results
write.csv(ego@result, "GO_Enrichment.csv", row.names = FALSE)


# Save annotated DESeq2 results
write.csv(resDF, file = "DESeq2_Annotated_Results.csv", row.names = TRUE)


getwd()


write.csv(resDF, file = "C:/Users/Ganga Supriya/OneDrive/Documents/DESeq2_Annotated_Results.csv", row.names = TRUE)

save.image(file = "my_workspace.RData")


