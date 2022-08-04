#1. plot proportion of the genes interacted with HIV gene
condition <- c("untreat", "short", "long", "ec")
percentage <- c((47/121)*100, (60/148)*100, (66/176)*100, (31/104)*100)
df_IS_interacted_HIV_gene <- data.frame(condition, percentage)

order_condition <- c("untreat", "short", "long", "ec")
df_IS_interacted_HIV_gene$condition <- factor(df_IS_interacted_HIV_gene$condition, levels = order_condition)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/Barplot_IS_interacted_HIV_gene.svg", height = 4, width = 2.8)
ggplot(df_IS_interacted_HIV_gene, aes(x = condition, y = percentage, fill = condition))+geom_bar(stat = "identity", color = "black")+theme_bw()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, vjust=1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

#2. plot clustering heatmap of the frequency of IS interacted with HIV genes in each physiological condition
HIV_genes <- c("Gag-pol", "Vif", "Vpr", "Tat", "Vpu", "Rev", "Env", "Nef")
untreat <- c(14/37, 7/37, 5/37, 13/37, 5/37, 3/37, 11/37, 6/37)
short <- c(21/51, 6/51, 6/51, 24/51, 0/51, 5/51, 16/51, 14/51)
long <- c(25/52, 3/52, 7/52, 18/52, 0/52, 3/52, 18/52, 12/52)
ec <- c(13/23, 1/23, 5/23, 7/23, 1/23, 2/23, 4/23, 7/23)

df_IS_interacted_HIV_gene_heatmap <- data.frame(HIV_genes, untreat, short, long, ec)
rownames(df_IS_interacted_HIV_gene_heatmap) <- df_IS_interacted_HIV_gene_heatmap$HIV_genes
df_IS_interacted_HIV_gene_heatmap$HIV_genes <- NULL
df_IS_interacted_HIV_gene_heatmap.mx <- data.matrix(df_IS_interacted_HIV_gene_heatmap, rownames.force = T)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/df_IS_interacted_HIV_gene_heatmap.mx.svg", height = 2.8, width = 4)
Heatmap(df_IS_interacted_HIV_gene_heatmap.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.5), column_title = "clinical condition", row_title = "HIV genes", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))

#3. import list of IS interacted/ not interacted with HIV genes
untreat_inter_hiv <- read.table("/media/chen/LaCie/IGH_backup/evoPath/df/IS_untreat_cART_as_hiv.csv", stringsAsFactors = F) %>% dplyr::rename(gene_name = V1) %>% dplyr::mutate(condition = "untreat", HIV_gene = "positive")
untreat_NOT_inter_hiv <- read.table("/media/chen/LaCie/IGH_backup/evoPath/df/IS_untreat_cART_nas_hiv.csv", stringsAsFactors = F) %>% dplyr::rename(gene_name = V1) %>% dplyr::mutate(condition = "untreat", HIV_gene = "negative")

short_inter_hiv <- read.table("/media/chen/LaCie/IGH_backup/evoPath/df/IS_short_cART_as_hiv.csv", stringsAsFactors = F) %>% dplyr::rename(gene_name = V1) %>% dplyr::mutate(condition = "short", HIV_gene = "positive")
short_NOT_inter_hiv <- read.table("/media/chen/LaCie/IGH_backup/evoPath/df/IS_short_cART_nas_hiv.csv", stringsAsFactors = F) %>% dplyr::rename(gene_name = V1) %>% dplyr::mutate(condition = "short", HIV_gene = "negative")

long_inter_hiv <- read.table("/media/chen/LaCie/IGH_backup/evoPath/df/IS_long_cART_as_hiv.csv", stringsAsFactors = F) %>% dplyr::rename(gene_name = V1) %>% dplyr::mutate(condition = "long", HIV_gene = "positive")
long_NOT_inter_hiv <- read.table("/media/chen/LaCie/IGH_backup/evoPath/df/IS_long_cART_nas_hiv.csv", stringsAsFactors = F) %>% dplyr::rename(gene_name = V1) %>% dplyr::mutate(condition = "long", HIV_gene = "negative")

tout_gene_list_inter_hiv <- dplyr::bind_rows(untreat_inter_hiv, short_inter_hiv, long_inter_hiv) %>% dplyr::select(gene_name) %>% unique()
tout_gene_list_NOT_inter_hiv <- dplyr::bind_rows(untreat_NOT_inter_hiv, short_NOT_inter_hiv, long_NOT_inter_hiv) %>% dplyr::select(gene_name) %>% unique()

#4. C7 enrichment analysis::tous
m_df <- msigdbr(species = "Homo sapiens")
m_tC7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% dplyr::select(gs_name, entrez_gene)

tout_gene_list_inter_hiv$gene_name <- as.character(tout_gene_list_inter_hiv$gene_name)
tout_gene_list_inter_hiv.enID <- bitr(tout_gene_list_inter_hiv$gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_tout_gene_list_inter_hiv <- data.frame(enricher(tout_gene_list_inter_hiv.enID$ENTREZID, TERM2GENE = m_tC7))


tout_gene_list_NOT_inter_hiv$gene_name <- as.character(tout_gene_list_NOT_inter_hiv$gene_name)
tout_gene_list_NOT_inter_hiv.enID <- bitr(tout_gene_list_NOT_inter_hiv$gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
msigdb_tout_gene_list_NOT_inter_hiv <- data.frame(enricher(tout_gene_list_NOT_inter_hiv.enID$ENTREZID, TERM2GENE = m_tC7))

## Ratio rechnen
overlap_GSE_ID <- dplyr::inner_join(msigdb_tout_gene_list_inter_hiv, msigdb_tout_gene_list_NOT_inter_hiv, by = "ID") %>% dplyr::select(ID) %>% unique() 
uniq_id_msigdb_tout_gene_list_inter_hiv <- dplyr::anti_join(msigdb_tout_gene_list_inter_hiv, overlap_GSE_ID, by = "ID")
write.table(uniq_id_msigdb_tout_gene_list_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/uniq_id_msigdb_tout_gene_list_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)
uniq_id_msigdb_tout_gene_list_NOT_inter_hiv <- dplyr::anti_join(msigdb_tout_gene_list_NOT_inter_hiv, overlap_GSE_ID, by = "ID")
write.table(uniq_id_msigdb_tout_gene_list_NOT_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/uniq_id_msigdb_tout_gene_list_NOT_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)

msigdb_tout_gene_list_inter_hiv$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_tout_gene_list_inter_hiv$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_tout_gene_list_inter_hiv$GeneRatio, perl=T))
msigdb_tout_gene_list_inter_hiv$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_tout_gene_list_inter_hiv$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_tout_gene_list_inter_hiv$BgRatio, perl=T))
 
msigdb_tout_gene_list_inter_hiv <- dplyr::mutate(msigdb_tout_gene_list_inter_hiv, Ratio = GeneRatio/BgRatio)

msigdb_tout_gene_list_NOT_inter_hiv$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_tout_gene_list_NOT_inter_hiv$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_tout_gene_list_NOT_inter_hiv$GeneRatio, perl=T))
msigdb_tout_gene_list_NOT_inter_hiv$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_tout_gene_list_NOT_inter_hiv$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_tout_gene_list_NOT_inter_hiv$BgRatio, perl=T))

msigdb_tout_gene_list_NOT_inter_hiv <- dplyr::mutate(msigdb_tout_gene_list_NOT_inter_hiv, Ratio = GeneRatio/BgRatio)

## heatmap
tout_GSE_ID <- dplyr::bind_rows(msigdb_tout_gene_list_inter_hiv, msigdb_tout_gene_list_NOT_inter_hiv) %>% dplyr::select(ID) %>% unique()

msigdb_tout_gene_list_inter_hiv_present <- dplyr::inner_join(msigdb_tout_gene_list_inter_hiv, tout_GSE_ID, by = "ID") %>% dplyr::select(ID, Ratio) %>% dplyr::rename(Interaction = Ratio)
msigdb_tout_gene_list_inter_hiv_absent <- dplyr::anti_join( tout_GSE_ID, msigdb_tout_gene_list_inter_hiv, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Interaction = 0)
msigdb_tout_gene_list_inter_hiv_combined <- dplyr::bind_rows(msigdb_tout_gene_list_inter_hiv_present, msigdb_tout_gene_list_inter_hiv_absent) %>% dplyr::arrange(desc(ID))

msigdb_tout_gene_list_NOT_inter_hiv_present <- dplyr::inner_join(msigdb_tout_gene_list_NOT_inter_hiv, tout_GSE_ID, by = "ID") %>% dplyr::select(ID, Ratio) %>% dplyr::rename(NOinteraction = Ratio)
msigdb_tout_gene_list_NOT_inter_hiv_absent <- dplyr::anti_join( tout_GSE_ID, msigdb_tout_gene_list_NOT_inter_hiv, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(NOinteraction = 0)
msigdb_tout_gene_list_NOT_inter_hiv_combined <- dplyr::bind_rows(msigdb_tout_gene_list_NOT_inter_hiv_present, msigdb_tout_gene_list_NOT_inter_hiv_absent) %>% dplyr::arrange(desc(ID))

msigdb_tout_gene_list_tous_hiv_combined <- dplyr::bind_cols(msigdb_tout_gene_list_inter_hiv_combined, msigdb_tout_gene_list_NOT_inter_hiv_combined) %>% dplyr::select(ID...1, Interaction, NOinteraction)
rownames(msigdb_tout_gene_list_tous_hiv_combined) <- msigdb_tout_gene_list_tous_hiv_combined$ID...1
msigdb_tout_gene_list_tous_hiv_combined$ID...1 <- NULL

msigdb_tout_gene_list_tous_hiv_combined.mx <- data.matrix(msigdb_tout_gene_list_tous_hiv_combined, rownames.force = T)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/msigdb_tout_gene_list_tous_hiv_combined.mx.svg", height = 5, width = 4)
Heatmap(msigdb_tout_gene_list_tous_hiv_combined.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.01), column_title = "clinical condition", row_title = "C7 GSE terms", show_row_names = F, show_column_dend = FALSE, column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

#5. C7 enrichment analysis::breakdown list
# C7 enrichment analysis
Get_entreID <- function(df) {
  df$gene_name <- as.character(df$gene_name)
df.enID <- bitr(df$gene_name, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
df.enID.m_tC7 <- data.frame(enricher(df.enID$ENTREZID, TERM2GENE = m_tC7))

return(df.enID.m_tC7)
}

msigdb_untreat_inter_hiv <- Get_entreID(untreat_inter_hiv)
write.table(msigdb_untreat_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/msigdb_untreat_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)
msigdb_untreat_NOT_inter_hiv <- Get_entreID(untreat_NOT_inter_hiv)
write.table(msigdb_untreat_NOT_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/msigdb_untreat_NOT_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)

msigdb_short_inter_hiv <- Get_entreID(short_inter_hiv)
write.table(msigdb_short_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/msigdb_short_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)
msigdb_short_NOT_inter_hiv <- Get_entreID(short_NOT_inter_hiv)
write.table(msigdb_short_NOT_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/msigdb_short_NOT_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)

msigdb_long_inter_hiv <- Get_entreID(long_inter_hiv)
write.table(msigdb_long_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/msigdb_long_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)
msigdb_long_NOT_inter_hiv <- Get_entreID(long_NOT_inter_hiv)
write.table(msigdb_long_NOT_inter_hiv, file = "/media/chen/LaCie/IGH_backup/evoPath/df/msigdb_long_NOT_inter_hiv.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# Ratio rechnen
Ratio_rechnen <- function(df) {
  df$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", df$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", df$GeneRatio, perl=T))
  df$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", df$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", df$BgRatio, perl=T))
 
df.ratio <- dplyr::mutate(df, Ratio = GeneRatio/BgRatio)

return(df.ratio)
}

msigdb_untreat_inter_hiv.ratio <- Ratio_rechnen(msigdb_untreat_inter_hiv)
msigdb_untreat_NOT_inter_hiv.ratio <- Ratio_rechnen(msigdb_untreat_NOT_inter_hiv)
msigdb_short_inter_hiv.ratio <- Ratio_rechnen(msigdb_short_inter_hiv)
msigdb_short_NOT_inter_hiv.ratio <- Ratio_rechnen(msigdb_short_NOT_inter_hiv)
msigdb_long_inter_hiv.ratio <- Ratio_rechnen(msigdb_long_inter_hiv)
msigdb_long_NOT_inter_hiv.ratio <- Ratio_rechnen(msigdb_long_NOT_inter_hiv)

# Heatmap machen
tout_GSE_ID_breakdown <- dplyr::bind_rows(msigdb_untreat_inter_hiv.ratio, msigdb_untreat_NOT_inter_hiv.ratio, msigdb_short_inter_hiv.ratio, msigdb_short_NOT_inter_hiv.ratio, msigdb_long_inter_hiv.ratio, msigdb_long_NOT_inter_hiv.ratio) %>% dplyr::select(ID) %>% unique()

Tous_GSE_ID_per_condition <- function(df, tous_GSE_ID) {
  df.present <- dplyr::inner_join(df, tous_GSE_ID, by = "ID") %>% dplyr::select(ID, Ratio)
  df.absent <- dplyr::anti_join( tous_GSE_ID, df, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
  df.combined <- dplyr::bind_rows(df.present, df.absent) %>% dplyr::arrange(desc(ID))
  
  return(df.combined)
}

msigdb_untreat_inter_hiv.ratio.tous_GSE <- Tous_GSE_ID_per_condition(msigdb_untreat_inter_hiv.ratio, tout_GSE_ID_breakdown) %>% dplyr::rename(Untreat_inter = Ratio)
msigdb_untreat_NOT_inter_hiv.ratio.tous_GSE <- Tous_GSE_ID_per_condition(msigdb_untreat_NOT_inter_hiv.ratio, tout_GSE_ID_breakdown) %>% dplyr::rename(Untreat_NOT_inter = Ratio)
msigdb_short_inter_hiv.ratio.tous_GSE <- Tous_GSE_ID_per_condition(msigdb_short_inter_hiv.ratio, tout_GSE_ID_breakdown) %>% dplyr::rename(Short_inter = Ratio)
msigdb_short_NOT_inter_hiv.ratio.tous_GSE <- Tous_GSE_ID_per_condition(msigdb_short_NOT_inter_hiv.ratio, tout_GSE_ID_breakdown) %>% dplyr::rename(Short_NOT_inter = Ratio)
msigdb_long_inter_hiv.ratio.tous_GSE <- Tous_GSE_ID_per_condition(msigdb_long_inter_hiv.ratio, tout_GSE_ID_breakdown) %>% dplyr::rename(Long_inter = Ratio)
msigdb_long_NOT_inter_hiv.ratio.tous_GSE <- Tous_GSE_ID_per_condition(msigdb_long_NOT_inter_hiv.ratio, tout_GSE_ID_breakdown) %>% dplyr::rename(Long_NOT_inter = Ratio)

Tous_msigdb_breakdown.ratio.tous_GSE <- dplyr::bind_cols(msigdb_untreat_inter_hiv.ratio.tous_GSE, msigdb_untreat_NOT_inter_hiv.ratio.tous_GSE, msigdb_short_inter_hiv.ratio.tous_GSE, msigdb_short_NOT_inter_hiv.ratio.tous_GSE, msigdb_long_inter_hiv.ratio.tous_GSE, msigdb_long_NOT_inter_hiv.ratio.tous_GSE) %>% dplyr::select(ID...1, Untreat_inter, Untreat_NOT_inter, Short_inter, Short_NOT_inter, Long_inter, Long_NOT_inter)

rownames(Tous_msigdb_breakdown.ratio.tous_GSE) <- Tous_msigdb_breakdown.ratio.tous_GSE$ID...1
Tous_msigdb_breakdown.ratio.tous_GSE$ID...1 <- NULL
Tous_msigdb_breakdown.ratio.tous_GSE.mx <- data.matrix(Tous_msigdb_breakdown.ratio.tous_GSE, rownames.force = T)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/Tous_msigdb_breakdown.ratio.tous_GSE.mx.svg", height = 5, width = 4)
Heatmap(Tous_msigdb_breakdown.ratio.tous_GSE.mx, name = "Frequency", rect_gp = gpar(col = "white", lwd = 0.1), column_title = "clinical condition", row_title = "C7 GSE terms", show_row_names = F, column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

#6. clustring heatmap cell types breakdown list
msigdb_cell_type_breakdown <- read.table("/media/chen/LaCie/MacBook_backup/Boulot/Projets/Moi/evoPATH/df/HIV_genes_interacted/MSigDB_cell_type_breakdown.tsv", header = T, stringsAsFactors = F)

# tous
msigdb_cell_type_breakdown_tous <- msigdb_cell_type_breakdown
rownames(msigdb_cell_type_breakdown_tous) <- msigdb_cell_type_breakdown_tous$Cell_type
msigdb_cell_type_breakdown_tous$Cell_type <- NULL
msigdb_cell_type_breakdown_tous.mx <- data.matrix(msigdb_cell_type_breakdown_tous, rownames.force = T)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/msigdb_cell_type_breakdown_tous.mx.svg", height = 4, width = 6)
Heatmap(msigdb_cell_type_breakdown_tous.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "clinical condition", row_title = "Cell types", show_row_names = T, column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

# breakdown
msigdb_cell_type_breakdown_nur <- msigdb_cell_type_breakdown[, c(1, 5:10)]
rownames(msigdb_cell_type_breakdown_nur) <- msigdb_cell_type_breakdown_nur$Cell_type
msigdb_cell_type_breakdown_nur$Cell_type <- NULL
msigdb_cell_type_breakdown_nur.mx <- data.matrix(msigdb_cell_type_breakdown_nur, rownames.force = T)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/msigdb_cell_type_breakdown_nur.mx.svg", height = 4, width = 5.5)
Heatmap(msigdb_cell_type_breakdown_nur.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "clinical condition", row_title = "Cell types", show_row_names = T, column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

#7. stacked plot breakdown
Get_df_for_stacked_plot <- function(df) {
  df.1 <- df %>% dplyr::select(Cell_type, ART_untreat) %>% dplyr::rename(percentage = ART_untreat) %>% dplyr::mutate(condition = "ART (ut)")
  df.2 <- df %>% dplyr::select(Cell_type, ART_short) %>% dplyr::rename(percentage = ART_short) %>% dplyr::mutate(condition = "ART (st)")
  df.3 <- df %>% dplyr::select(Cell_type, ART_long) %>% dplyr::rename(percentage = ART_long) %>% dplyr::mutate(condition = "ART (lt)")
  df.4 <- df %>% dplyr::select(Cell_type, ART_untreat_hiv) %>% dplyr::rename(percentage = ART_untreat_hiv) %>% dplyr::mutate(condition = "ART (ut; +HIV)")
  df.5 <- df %>% dplyr::select(Cell_type, ART_untreat_NOhiv) %>% dplyr::rename(percentage = ART_untreat_NOhiv) %>% dplyr::mutate(condition = "ART (ut; -HIV)")
  df.6 <- df %>% dplyr::select(Cell_type, ART_short_hiv) %>% dplyr::rename(percentage = ART_short_hiv) %>% dplyr::mutate(condition = "ART (st; +HIV)")
  df.7 <- df %>% dplyr::select(Cell_type, ART_short_NOhiv) %>% dplyr::rename(percentage = ART_short_NOhiv) %>% dplyr::mutate(condition = "ART (st; -HIV)")
  df.8 <- df %>% dplyr::select(Cell_type, ART_long_hiv) %>% dplyr::rename(percentage = ART_long_hiv) %>% dplyr::mutate(condition = "ART (lt; +HIV)")
  df.9 <- df %>% dplyr::select(Cell_type, ART_long_NOhiv) %>% dplyr::rename(percentage = ART_long_NOhiv) %>% dplyr::mutate(condition = "ART (lt; -HIV)")
  
  df.tous <- dplyr::bind_rows(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8, df.9)
  
  deposer_cond <- c("ART (ut)", "ART (ut; +HIV)", "ART (ut; -HIV)", "ART (st)", "ART (st; +HIV)", "ART (st; -HIV)", "ART (lt)", "ART (lt; +HIV)", "ART (lt; -HIV)")  
  
  df.tous$condition <- factor(df.tous$condition, levels = deposer_cond)
  
  return(df.tous)
}

msigdb_cell_type_breakdown_stacked_bar_df <- Get_df_for_stacked_plot(msigdb_cell_type_breakdown)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/msigdb_cell_type_breakdown_stacked_bar_df.svg", height = 6, width = 6)  
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)
ggplot(msigdb_cell_type_breakdown_stacked_bar_df, aes(x = condition, y = percentage, fill = Cell_type))+geom_bar(stat = "identity", position = "stack")+scale_fill_manual(values = c(mycolors))+theme_bw()+ylim(0, 1.5)+labs(fill = "Cell Type")+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, vjust=1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

### ***chi-sequre test*** ###
chisq.test(msigdb_cell_type_breakdown[, 2], msigdb_cell_type_breakdown[, 5])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 2] and msigdb_cell_type_breakdown[, 5]
#X-squared = 17.231, df = 4, p-value = 0.001743
**************************************************************************************************
chisq.test(msigdb_cell_type_breakdown[, 3], msigdb_cell_type_breakdown[, 7])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 3] and msigdb_cell_type_breakdown[, 7]
#X-squared = 81.659, df = 48, p-value = 0.001748
**************************************************************************************************
chisq.test(msigdb_cell_type_breakdown[, 7], msigdb_cell_type_breakdown[, 8])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 7] and msigdb_cell_type_breakdown[, 8]
#X-squared = 41.697, df = 24, p-value = 0.01395
**************************************************************************************************
chisq.test(msigdb_cell_type_breakdown[, 4], msigdb_cell_type_breakdown[, 9])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 4] and msigdb_cell_type_breakdown[, 9]
#X-squared = 94.476, df = 70, p-value = 0.0273
**************************************************************************************************
chisq.test(msigdb_cell_type_breakdown[, 9], msigdb_cell_type_breakdown[, 10])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 9] and msigdb_cell_type_breakdown[, 10]
#X-squared = 32, df = 14, p-value = 0.004006
**************************************************************************************************  
chisq.test(msigdb_cell_type_breakdown[, 3], msigdb_cell_type_breakdown[, 8])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 3] and msigdb_cell_type_breakdown[, 8]
#X-squared = 43.273, df = 18, p-value = 0.0007322
**************************************************************************************************    
chisq.test(msigdb_cell_type_breakdown[, 4], msigdb_cell_type_breakdown[, 10])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 4] and msigdb_cell_type_breakdown[, 10]
#X-squared = 32, df = 20, p-value = 0.0433
**************************************************************************************************    
chisq.test(msigdb_cell_type_breakdown[, 7], msigdb_cell_type_breakdown[, 9])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 7] and msigdb_cell_type_breakdown[, 9]
#X-squared = 81.371, df = 56, p-value = 0.01502
**************************************************************************************************
chisq.test(msigdb_cell_type_breakdown[, 5], msigdb_cell_type_breakdown[, 9])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 5] and msigdb_cell_type_breakdown[, 9]
#X-squared = 15.473, df = 14, p-value = 0.3466
**************************************************************************************************
chisq.test(msigdb_cell_type_breakdown[, 8], msigdb_cell_type_breakdown[, 10])

**************************************************************************************************
#Pearson's Chi-squared test

#data:  msigdb_cell_type_breakdown[, 8] and msigdb_cell_type_breakdown[, 10]
#X-squared = 8.2797, df = 6, p-value = 0.2183
**************************************************************************************************  

#8. r cluster heatmap cytokine breakdown list
msigdb_cytokine_breakdown <- read.table("/media/chen/LaCie/MacBook_backup/Boulot/Projets/Moi/evoPATH/df/HIV_genes_interacted/MSigDB_cytokines_breakdown.tsv", header = T, stringsAsFactors = F)

row.names(msigdb_cytokine_breakdown) <- msigdb_cytokine_breakdown$Cytokine
msigdb_cytokine_breakdown$Cytokine <- NULL
msigdb_cytokine_breakdown.mx <- data.matrix(msigdb_cytokine_breakdown, rownames.force = T)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/msigdb_cytokine_breakdown.mx.svg", height = 5, width = 4.8)
Heatmap(msigdb_cytokine_breakdown.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "clinical condition", row_title = "Cytokines", show_row_names = T, column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()
