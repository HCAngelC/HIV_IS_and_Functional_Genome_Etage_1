#1. Run analysis
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

#2. selection key columns & compute ratio
msigdb_c7_cART_short.key <- msigdb_c7_cART_short %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART short")
msigdb_c7_cART_short.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$GeneRatio, perl=T))
msigdb_c7_cART_short.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$BgRatio, perl=T))

msigdb_c7_cART_long.key <- msigdb_c7_cART_long %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART long")
msigdb_c7_cART_long.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$GeneRatio, perl=T))
msigdb_c7_cART_long.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$BgRatio, perl=T))

msigdb_c7_cART_short.key <- msigdb_c7_cART_short.key %>% dplyr::mutate(genes = "all")
msigdb_c7_cART_short_sans_cancer_specific_type.key <- msigdb_c7_cART_short_sans_cancer_specific_type %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART short", genes = "removal cancer specific type")
msigdb_c7_cART_short_sans_cancer_specific_type.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_cancer_specific_type.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_cancer_specific_type.key$GeneRatio, perl=T))
msigdb_c7_cART_short_sans_cancer_specific_type.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_cancer_specific_type.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_cancer_specific_type.key$BgRatio, perl=T))

msigdb_c7_cART_short_sans_signal_transduction.key <- msigdb_c7_cART_short_sans_signal_transduction %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART short", genes = "removal signal transduction")
msigdb_c7_cART_short_sans_signal_transduction.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_signal_transduction.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_signal_transduction.key$GeneRatio, perl=T))
msigdb_c7_cART_short_sans_signal_transduction.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_signal_transduction.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_signal_transduction.key$BgRatio, perl=T))

msigdb_c7_cART_short_sans_immune_system.key <- msigdb_c7_cART_short_sans_immune_system %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART short", genes = "removal immune system")
msigdb_c7_cART_short_sans_immune_system.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_immune_system.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_immune_system.key$GeneRatio, perl=T))
msigdb_c7_cART_short_sans_immune_system.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_immune_system.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_immune_system.key$BgRatio, perl=T))

msigdb_c7_cART_short_sans_infectious_disease_viral.key <- msigdb_c7_cART_short_sans_infectious_disease_viral %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART short", genes = "removal infectious disease viral")
msigdb_c7_cART_short_sans_infectious_disease_viral.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_infectious_disease_viral.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_infectious_disease_viral.key$GeneRatio, perl=T))
msigdb_c7_cART_short_sans_infectious_disease_viral.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short_sans_infectious_disease_viral.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short_sans_infectious_disease_viral.key$BgRatio, perl=T))

msigdb_c7_cART_long.key <- msigdb_c7_cART_long.key %>% dplyr::mutate(genes = "all")
msigdb_c7_cART_long_sans_cancer_specific_type.key <- msigdb_c7_cART_long_sans_cancer_specific_type %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART long", genes = "removal cancer specific type")
msigdb_c7_cART_long_sans_cancer_specific_type.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_cancer_specific_type.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_cancer_specific_type.key$GeneRatio, perl=T))
msigdb_c7_cART_long_sans_cancer_specific_type.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_cancer_specific_type.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_cancer_specific_type.key$BgRatio, perl=T))

msigdb_c7_cART_long_sans_signal_transduction.key <- msigdb_c7_cART_long_sans_signal_transduction %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART long", genes = "removal signal transduction")
msigdb_c7_cART_long_sans_signal_transduction.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_signal_transduction.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_signal_transduction.key$GeneRatio, perl=T))
msigdb_c7_cART_long_sans_signal_transduction.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_signal_transduction.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_signal_transduction.key$BgRatio, perl=T))

msigdb_c7_cART_long_sans_immune_system.key <- msigdb_c7_cART_long_sans_immune_system %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART long", genes = "removal immune system")
msigdb_c7_cART_long_sans_immune_system.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_immune_system.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_immune_system.key$GeneRatio, perl=T))
msigdb_c7_cART_long_sans_immune_system.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_immune_system.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_immune_system.key$BgRatio, perl=T))

msigdb_c7_cART_long_sans_infectious_disease_viral.key <- msigdb_c7_cART_long_sans_infectious_disease_viral %>% dplyr::select(ID, GeneRatio, BgRatio, p.adjust) %>% dplyr::mutate(cond = "cART long", genes = "removal infectious disease viral")
msigdb_c7_cART_long_sans_infectious_disease_viral.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_infectious_disease_viral.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_infectious_disease_viral.key$GeneRatio, perl=T))
msigdb_c7_cART_long_sans_infectious_disease_viral.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long_sans_infectious_disease_viral.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long_sans_infectious_disease_viral.key$BgRatio, perl=T))

df_cART_short_long_tous_et_removel_top4KEGG <- dplyr::bind_rows(msigdb_c7_cART_short.key, msigdb_c7_cART_short_sans_cancer_specific_type.key, msigdb_c7_cART_short_sans_signal_transduction.key, msigdb_c7_cART_short_sans_immune_system.key, msigdb_c7_cART_short_sans_infectious_disease_viral.key, msigdb_c7_cART_long.key, msigdb_c7_cART_long_sans_cancer_specific_type.key, msigdb_c7_cART_long_sans_signal_transduction.key, msigdb_c7_cART_long_sans_immune_system.key, msigdb_c7_cART_long_sans_infectious_disease_viral.key) 

deposer_cond = c("cART short", "cART long")
deposer_genes = c("all", "removal cancer specific type", "removal signal transduction", "removal immune system", "removal infectious disease viral")
comparaison_genes <- list(c("all", "removal cancer specific type"), c("all", "removal signal transduction"), c("all", "removal immune system"), c("all", "removal infectious disease viral"))

df_cART_short_long_tous_et_removel_top4KEGG$cond <- factor(df_cART_short_long_tous_et_removel_top4KEGG$cond, levels = deposer_cond)
df_cART_short_long_tous_et_removel_top4KEGG$genes <- factor(df_cART_short_long_tous_et_removel_top4KEGG$genes, levels = deposer_genes)

#3. plot p.adj
pdf("/media/chen/DATA/evoPath/Abb/df_cART_short_long_tous_et_removel_top4KEGG_p.adjust.pdf", height = 5.2, width = 5.2) 
ggplot(df_cART_short_long_tous_et_removel_top4KEGG, aes(x = genes, y = p.adjust, fill = genes))+geom_violin(trim = F)+geom_boxplot(width = 0.1, fill = "white")+scale_fill_brewer(palette = "Blues")+stat_compare_means(comparisons = comparaison_genes, label = "p.signif")+facet_grid(. ~ cond)+theme_bw()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()
