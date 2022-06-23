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

#4. plot heatmap
df_cytokines_removal_top4_KEGG <- read.table("/cytokines_removal_top4_KEGG.csv", header = T, stringsAsFactors = F)
row.names(df_cytokines_removal_top4_KEGG) <- df_cytokines_removal_top4_KEGG$Cytokines
df_cytokines_removal_top4_KEGG$Cytokines <- NULL
df_cytokines_removal_top4_KEGG.mx <- data.matrix(df_cytokines_removal_top4_KEGG, rownames.force = T)

Heatmap(df_cytokines_removal_top4_KEGG.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.5), column_title = "clinical condition", row_title = "Cytokines", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))

#4. plot Correlation plot
# Obtenir le triangle inférieur
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Obtenir le triangle supérieur
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

# short cART
df_cytokines_removal_top4_KEGG.mx.short.cor <- cor(df_cytokines_removal_top4_KEGG.mx[,1:5], method = c("pearson"))
upper_df_cytokines_removal_top4_KEGG.mx.short.cor <- get_upper_tri(df_cytokines_removal_top4_KEGG.mx.short.cor)
melted_upper_df_cytokines_removal_top4_KEGG.mx.short.cor <- melt(upper_df_cytokines_removal_top4_KEGG.mx.short.cor, na.rm = T)

ggplot(data = melted_upper_df_cytokines_removal_top4_KEGG.mx.short.cor, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, limit = c(0, 1), space = "Lab", name = "Pearson\nCorrelation")+theme_minimal()+theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+coord_fixed()

# output
                       Var1                     Var2     value
1            cART_short_all           cART_short_all 1.0000000
6            cART_short_all     cART_short_rm_cancer 0.7888106
7      cART_short_rm_cancer     cART_short_rm_cancer 1.0000000
11           cART_short_all      cART_short_rm_immue 0.6179144
12     cART_short_rm_cancer      cART_short_rm_immue 0.2901294
13      cART_short_rm_immue      cART_short_rm_immue 1.0000000
16           cART_short_all     cART_short_rm_signal 0.7888106
17     cART_short_rm_cancer     cART_short_rm_signal 1.0000000
18      cART_short_rm_immue     cART_short_rm_signal 0.2901294
19     cART_short_rm_signal     cART_short_rm_signal 1.0000000
21           cART_short_all cART_short_rm_infectious 0.8874120
22     cART_short_rm_cancer cART_short_rm_infectious 0.8888889
23      cART_short_rm_immue cART_short_rm_infectious 0.4497006
24     cART_short_rm_signal cART_short_rm_infectious 0.8888889
25 cART_short_rm_infectious cART_short_rm_infectious 1.0000000

# long cART
df_cytokines_removal_top4_KEGG.mx.long.cor <- cor(df_cytokines_removal_top4_KEGG.mx[,6:10], method = c("pearson"))
upper_df_cytokines_removal_top4_KEGG.mx.long.cor <- get_upper_tri(df_cytokines_removal_top4_KEGG.mx.long.cor)
melted_upper_df_cytokines_removal_top4_KEGG.mx.long.cor <- melt(upper_df_cytokines_removal_top4_KEGG.mx.long.cor, na.rm = T)

ggplot(data = melted_upper_df_cytokines_removal_top4_KEGG.mx.long.cor, aes(Var2, Var1, fill = value))+geom_tile(color = "white")+scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, limit = c(0, 1), space = "Lab", name = "Pearson\nCorrelation")+theme_minimal()+theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))+coord_fixed()

# output
                      Var1                    Var2     value
1            cART_long_all           cART_long_all 1.0000000
6            cART_long_all     cART_long_rm_cancer 0.9308937
7      cART_long_rm_cancer     cART_long_rm_cancer 1.0000000
11           cART_long_all      cART_long_rm_immue 0.8082827
12     cART_long_rm_cancer      cART_long_rm_immue 0.7634454
13      cART_long_rm_immue      cART_long_rm_immue 1.0000000
16           cART_long_all     cART_long_rm_signal 0.6466800
17     cART_long_rm_cancer     cART_long_rm_signal 0.6598226
18      cART_long_rm_immue     cART_long_rm_signal 0.8215982
19     cART_long_rm_signal     cART_long_rm_signal 1.0000000
21           cART_long_all cART_long_rm_infectious 0.7160235
22     cART_long_rm_cancer cART_long_rm_infectious 0.7885199
23      cART_long_rm_immue cART_long_rm_infectious 0.7546869
24     cART_long_rm_signal cART_long_rm_infectious 0.6620612
25 cART_long_rm_infectious cART_long_rm_infectious 1.0000000
