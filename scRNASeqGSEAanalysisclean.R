library(Seurat)
library(presto)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)

# identify idents in merged object
table(mergefin$orig.ident)

# find markers
CMU.markers823_2 <- FindMarkers(mergefin, group.by = "orig.ident",only.pos = FALSE, ident.1 = "CMU823", ident.2 = NULL, test.use = "DESeq2", max.cells.per.ident = 50)
original_gene_list4CMU <- CMU.markers823_2$avg_log2FC
names(original_gene_list4CMU) <- row.names(CMU.markers823_2)
gene_list4CMU <- na.omit(original_gene_list4CMU)
gene_list4CMU <- sort(original_gene_list4CMU, decreasing = TRUE)
gseCMU1 <- gseGO(geneList=gene_list4CMU,
+                  ont ="ALL",
+                  keyType = "SYMBOL",
+                  nPerm = 10000,
+                  minGSSize = 3,
+                  maxGSSize = 800,
+                  pvalueCutoff = 0.05,
+                  verbose = TRUE,
+                  OrgDb = huorganism,
+                  pAdjustMethod = "none")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)