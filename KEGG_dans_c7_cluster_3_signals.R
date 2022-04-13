# 1. Import df list genes cluster 3
c7_long_list_gene_TNF <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_TNF.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_TNF.id <- bitr(c7_long_list_gene_TNF$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_CXCL4 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_CXCL4.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_CXCL4.id <- bitr(c7_long_list_gene_CXCL4$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IL18 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IL18.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IL18.id <- bitr(c7_long_list_gene_IL18$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IL15 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IL15.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IL15.id <- bitr(c7_long_list_gene_IL15$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IL7 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IL7.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IL7.id <- bitr(c7_long_list_gene_IL7$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IL1 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IL1.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IL1.id <- bitr(c7_long_list_gene_IL1$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IL2 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IL2.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IL2.id <- bitr(c7_long_list_gene_IL2$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IFNA <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IFNA.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IFNA.id <- bitr(c7_long_list_gene_IFNA$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

c7_long_list_gene_IL6 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/c7_list_genes/c7_long_IL6.rtf", header = F, stringsAsFactors = T) %>% unique()
c7_long_list_gene_IL6.id <- bitr(c7_long_list_gene_IL6$V1, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL"), OrgDb="org.Hs.eg.db")

# 2. KEGG pathway over-representation analysis on list genes cluster 3
kegg_c7_long_list_gene_TNF <- data.frame(enrichKEGG(gene = c7_long_list_gene_TNF.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5))
write.table(kegg_c7_long_list_gene_TNF, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_TNF.txt", row.names = F, col.names = T, sep = "\t", quote = F) ## No pathway is enriched.

kegg_c7_long_list_gene_CXCL4 <- data.frame(enrichKEGG(gene = c7_long_list_gene_CXCL4.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_CXCL4, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_CXCL4.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IL18 <- data.frame(enrichKEGG(gene = c7_long_list_gene_IL18.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IL18, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IL18.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IL15 <- data.frame(enrichKEGG(gene = c7_long_list_gene_IL15.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IL15, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IL15.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IL7 <- data.frame(enrichKEGG(gene = c7_long_list_gene_IL7.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IL7, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IL7.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IL1 <- data.frame(enrichKEGG(gene = c7_long_list_gene_IL1.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IL1, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IL1.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IL2 <- data.frame(enrichKEGG(gene = c7_long_list_gene_IL2.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IL2, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IL2.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IFNA <- data.frame(enrichKEGG(gene = c7_long_list_gene_IFNA.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IFNA, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IFNA.txt", row.names = F, col.names = T, sep = "\t", quote = F)

kegg_c7_long_list_gene_IL6 <- data.frame(enrichKEGG(gene = c7_long_list_gene_IL6.id$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)) 
write.table(kegg_c7_long_list_gene_IL6, file = "/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_list_gene_IL6.txt", row.names = F, col.names = T, sep = "\t", quote = F) ## No pathway is enriched.
