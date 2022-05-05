### ***Heapmap (tous des conditions)*** ###
df_c7_tous_cell_type <- read.table("/media/chen/DATA/evoPath/df/C7/C7_msigdb_cell_types.tsv", header = T, stringsAsFactors = F)

df_c7_tous_cell_type.v <- df_c7_tous_cell_type
row.names(df_c7_tous_cell_type.v) <- df_c7_tous_cell_type.v$Cell_type
df_c7_tous_cell_type.v$Cell_type <- NULL
df_c7_tous_cell_type.v.mx <- data.matrix(df_c7_tous_cell_type.v, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/df_c7_tous_cell_type.v.mx.pdf", height = 5, width = 5.5)  
Heatmap(df_c7_tous_cell_type.v.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Clinical condition", row_title = "Cell type", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/df_c7_tous_cell_type.v.mx.svg", height = 5, width = 5.5)  
Heatmap(df_c7_tous_cell_type.v.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Clinical condition", row_title = "Cell type", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

### ***stacked barplot (sans EC)*** ###
df_c7_tous_cell_type.cART_untreat <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, cART_untreat) %>% dplyr::rename(percentage = cART_untreat) %>% dplyr::mutate(condition = "cART (ut)")
df_c7_tous_cell_type.cART_short <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, cART_short) %>% dplyr::rename(percentage = cART_short) %>% dplyr::mutate(condition = "cART (st)")
df_c7_tous_cell_type.cART_long <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, cART_long) %>% dplyr::rename(percentage = cART_long) %>% dplyr::mutate(condition = "cART (lt)")

df_c7_tous_cell_type_sansEC_barplot <- dplyr::bind_rows(df_c7_tous_cell_type.cART_untreat, df_c7_tous_cell_type.cART_short, df_c7_tous_cell_type.cART_long)

deposer_cond_sans_EC <- c("cART (ut)", "cART (st)", "cART (lt)")
df_c7_tous_cell_type_sansEC_barplot$condition <- factor(df_c7_tous_cell_type_sansEC_barplot$condition, levels = deposer_cond_sans_EC)

pdf("/media/chen/DATA/evoPath/Abb/df_c7_tous_cell_type_sans_EC_barplot.pdf", height = 9, width = 4.2)  
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)
ggplot(df_c7_tous_cell_type_sansEC_barplot, aes(x = condition, y = percentage, fill = Cell_type))+geom_bar(stat = "identity", position = "stack")+scale_fill_manual(values = addalpha(mycolors))+theme_bw()+scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1), breaks = seq(0,1, 0.1))+labs(fill = "Cell Type")+geom_line(aes(color = Cell_type, group = Cell_type, y = percentage), size = 1)+geom_point()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, vjust=1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

### ***chi-sequre test*** ###
chisq.test(df_c7_tous_cell_type.cART_untreat$percentage, df_c7_tous_cell_type.cART_short$percentage)

**************************************************************************************************
#Pearson's Chi-squared test

#data:  df_c7_tous_cell_type.cART_untreat$percentage and df_c7_tous_cell_type.cART_short$percentage
#X-squared = 22.4, df = 12, p-value = 0.03327
**************************************************************************************************

chisq.test(df_c7_tous_cell_type.cART_short$percentage, df_c7_tous_cell_type.cART_long$percentage)

**************************************************************************************************
#Pearson's Chi-squared test

#data:  df_c7_tous_cell_type.cART_short$percentage and df_c7_tous_cell_type.cART_long$percentage
#X-squared = 89.481, df = 60, p-value = 0.008108
**************************************************************************************************

chisq.test(df_c7_tous_cell_type.cART_untreat$percentage, df_c7_tous_cell_type.cART_long$percentage)

**************************************************************************************************
#Pearson's Chi-squared test

#data:  df_c7_tous_cell_type.cART_untreat$percentage and df_c7_tous_cell_type.cART_long$percentage
#X-squared = 28.8, df = 20, p-value = 0.09177
**************************************************************************************************

### ***stacked barplot (tous des scenarios)*** ###
df_c7_tous_cell_type.cART_untreat <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, cART_untreat) %>% dplyr::rename(percentage = cART_untreat) %>% dplyr::mutate(condition = "cART (ut)")
df_c7_tous_cell_type.cART_short <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, cART_short) %>% dplyr::rename(percentage = cART_short) %>% dplyr::mutate(condition = "cART (st)")
df_c7_tous_cell_type.cART_long <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, cART_long) %>% dplyr::rename(percentage = cART_long) %>% dplyr::mutate(condition = "cART (lt)")
df_c7_tous_cell_type.EC_long <- df_c7_tous_cell_type %>% dplyr::select(Cell_type, EC_long) %>% dplyr::rename(percentage = EC_long) %>% dplyr::mutate(condition = "Long-term EC")

df_c7_tous_cell_type_barplot <- dplyr::bind_rows(df_c7_tous_cell_type.cART_untreat, df_c7_tous_cell_type.cART_short, df_c7_tous_cell_type.cART_long, df_c7_tous_cell_type.EC_long)

deposer_cond <- c("cART (ut)", "cART (st)", "cART (lt)", "Long-term EC")
df_c7_tous_cell_type_barplot$condition <- factor(df_c7_tous_cell_type_barplot$condition, levels = deposer_cond)

pdf("/media/chen/DATA/evoPath/Abb/df_c7_tous_cell_type_barplot.pdf", height = 7, width = 4.8)  
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)
ggplot(df_c7_tous_cell_type_barplot, aes(x = condition, y = percentage, fill = Cell_type))+geom_bar(stat = "identity", position = "stack")+scale_fill_manual(values = addalpha(mycolors))+theme_bw()+scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1), breaks = seq(0,1, 0.1))+labs(fill = "Cell Type")+geom_line(aes(color = Cell_type, group = Cell_type, y = percentage), size = 1)+geom_point()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, vjust=1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()
