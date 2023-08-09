library(Seurat)
library(presto)
library(devtools)
install_github('immunogenomics/presto')
library(clusterProfiler)
library(enrichplot)
library(ggplot2)


BiocManager::install(huorganism, character.only = TRUE)
ids4PP2 <-bitr(names(original_gene_list4PP), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=huorganism)
dedup_ids4PP2 = ids[!duplicated(ids4PP2[c("ENSEMBL")]),]
dedup_ids4PP2 = ids4PP2[!duplicated(ids4PP2[c("ENSEMBL")]),]
dedup_ids4PP2 = ids4PP2[!duplicated(ids4PP2[c("SYMBOL")]),]
dotplot(kk2PP, showCategory = 10, title = "Enriched Pathways", split = ".sign") + facet_grid(.~.sign)
library(pathview)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="dme000010", species = kegg_organism)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="000010", species = kegg_organism)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="01100", species = kegg_organism)
pathview(gene.data = kegg_gene_list, pathway.id = "hsa01100", species = kegg_organism, limit = list(gene=5, cpd=1))
pathview(gene.data = kegg_gene_list, pathway.id = "hsa00010", species = kegg_organism, limit = list(gene=5, cpd=1))
pathview(gene.data = kegg_gene_list, pathway.id = "hsa00190", species = kegg_organism, limit = list(gene=5, cpd=1))
pathview(gene.data = kegg_gene_list, pathway.id = "hsa00020", species = kegg_organism, limit = list(gene=5, cpd=1))