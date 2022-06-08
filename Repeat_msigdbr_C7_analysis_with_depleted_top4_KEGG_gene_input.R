# sans cancer specific type
cART_short_sans_cancer_specific_type <- dplyr::anti_join(cART_short, list_gene_cancer_specific_type, by = "Gene_name")
cART_short_sans_cancer_specific_type$Gene_name <- as.character(cART_short_sans_cancer_specific_type$Gene_name)
cART_short_sans_cancer_specific_type.enID <- bitr(cART_short_sans_cancer_specific_type$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_short_sans_cancer_specific_type <- data.frame(enricher(cART_short_sans_cancer_specific_type.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_short_sans_cancer_specific_type, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_short_sans_cancer_specific_type.txt", row.names = F, col.names = T, sep = "\t", quote = F)

cART_long_sans_cancer_specific_type <- dplyr::anti_join(cART_long, list_gene_cancer_specific_type, by = "Gene_name")
cART_long_sans_cancer_specific_type$Gene_name <- as.character(cART_long_sans_cancer_specific_type$Gene_name)
cART_long_sans_cancer_specific_type.enID <- bitr(cART_long_sans_cancer_specific_type$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_long_sans_cancer_specific_type <- data.frame(enricher(cART_long_sans_cancer_specific_type.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_long_sans_cancer_specific_type, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_long_sans_cancer_specific_type.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# sans signal transduction
cART_short_sans_signal_transduction <- dplyr::anti_join(cART_short, list_gene_signal_transduction, by = "Gene_name")
cART_short_sans_signal_transduction$Gene_name <- as.character(cART_short_sans_signal_transduction$Gene_name)
cART_short_sans_signal_transduction.enID <- bitr(cART_short_sans_signal_transduction$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_short_sans_signal_transduction <- data.frame(enricher(cART_short_sans_signal_transduction.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_short_sans_signal_transduction, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_short_sans_signal_transduction.txt", row.names = F, col.names = T, sep = "\t", quote = F)

cART_long_sans_signal_transduction <- dplyr::anti_join(cART_long, list_gene_signal_transduction, by = "Gene_name")
cART_long_sans_signal_transduction$Gene_name <- as.character(cART_long_sans_signal_transduction$Gene_name)
cART_long_sans_signal_transduction.enID <- bitr(cART_long_sans_signal_transduction$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_long_sans_signal_transduction <- data.frame(enricher(cART_long_sans_signal_transduction.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_long_sans_signal_transduction, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_long_sans_signal_transduction.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# sans immune system
cART_short_sans_immune_system <- dplyr::anti_join(cART_short, list_gene_immune_system, by = "Gene_name")
cART_short_sans_immune_system$Gene_name <- as.character(cART_short_sans_immune_system$Gene_name)
cART_short_sans_immune_system.enID <- bitr(cART_short_sans_immune_system$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_short_sans_immune_system <- data.frame(enricher(cART_short_sans_immune_system.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_short_sans_immune_system, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_short_sans_immune_system.txt", row.names = F, col.names = T, sep = "\t", quote = F)

cART_long_sans_immune_system <- dplyr::anti_join(cART_long, list_gene_immune_system, by = "Gene_name")
cART_long_sans_immune_system$Gene_name <- as.character(cART_long_sans_immune_system$Gene_name)
cART_long_sans_immune_system.enID <- bitr(cART_long_sans_immune_system$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_long_sans_immune_system <- data.frame(enricher(cART_long_sans_immune_system.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_long_sans_immune_system, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_long_sans_immune_system.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# sans infectious disease viral
cART_short_sans_infectious_disease_viral <- dplyr::anti_join(cART_short, list_gene_infectious_disease_viral, by = "Gene_name")
cART_short_sans_infectious_disease_viral$Gene_name <- as.character(cART_short_sans_infectious_disease_viral$Gene_name)
cART_short_sans_infectious_disease_viral.enID <- bitr(cART_short_sans_infectious_disease_viral$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_short_sans_infectious_disease_viral <- data.frame(enricher(cART_short_sans_infectious_disease_viral.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_short_sans_infectious_disease_viral, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_short_sans_infectious_disease_viral.txt", row.names = F, col.names = T, sep = "\t", quote = F)

cART_long_sans_infectious_disease_viral <- dplyr::anti_join(cART_long, list_gene_infectious_disease_viral, by = "Gene_name")
cART_long_sans_infectious_disease_viral$Gene_name <- as.character(cART_long_sans_infectious_disease_viral$Gene_name)
cART_long_sans_infectious_disease_viral.enID <- bitr(cART_long_sans_infectious_disease_viral$Gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_c7_cART_long_sans_infectious_disease_viral <- data.frame(enricher(cART_long_sans_infectious_disease_viral.enID$ENTREZID, TERM2GENE = m_tC7))
write.table(msigdb_c7_cART_long_sans_infectious_disease_viral, file = "/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/msigdb_C7_cART_long_sans_infectious_disease_viral.txt", row.names = F, col.names = T, sep = "\t", quote = F)
