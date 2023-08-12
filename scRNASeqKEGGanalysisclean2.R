library(Seurat)
library(presto)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)

# Load and identify idents in merged object
load("Z:/Andre/mergefinprocessed1.RData")
table(mergefin$orig.ident)

# Find markers
CMU.markers823_2 <- FindMarkers(mergefin, group.by = "orig.ident",only.pos = FALSE, ident.1 = "CMU823", ident.2 = NULL, test.use = "DESeq2", max.cells.per.ident = 50)
original_gene_list4CMU <- CMU.markers823_2$avg_log2FC
names(original_gene_list4CMU) <- row.names(CMU.markers823_2)
gene_list4CMU <- na.omit(original_gene_list4CMU)
gene_list4CMU <- sort(original_gene_list4CMU, decreasing = TRUE)

# Transform gene IDs to Entrez IDs
ids4CMU1 <-bitr(names(original_gene_list4CMU), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=huorganism)

# Initialize list
my_list_6 <- list()

# Find indices for successfully mapped vector elements
newCMUindices <- match(ids4CMU1$SYMBOL, names(gene_list4CMU))

# Create new vector with mapped reads and IDs
for(i in 1:length(newCMUindices)) {
my_list_6[[i]] <- gene_list4CMU[[newCMUindices[[i]] ]]
}

# New vectors without Entrez ID transformation
newnumsCMU <- unlist(my_list_6)
CMU_KEGG_genelist <- newnumsCMU

# New vectors with Entrez ID transformation
ids4CMU1_1 <-bitr(names(gene_list4CMU), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=huorganism)
names(CMU_KEGG_genelist) <- ids4CMU1_1$ENTREZID

# View pathway analysis plot for desired pathway
pathview(gene.data = CMU_KEGG_genelist,
pathway.id = demo.paths$sel.paths[1],
species = "hsa",
out.suffix = "gse16873",
kegg.native = T)

