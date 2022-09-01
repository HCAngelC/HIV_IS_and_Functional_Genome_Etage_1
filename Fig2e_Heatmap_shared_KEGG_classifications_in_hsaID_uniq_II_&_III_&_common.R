## version 1 ##
group_II_uniq <- read.table("/media/chen/DATA/hasID_Cluster_II_uniq_count.agg.txt", header = T, stringsAsFactors = F) %>% dplyr::select(KEGG_subcat, count)
group_III_uniq <- read.table("/media/chen/DATA/hasID_Cluster_III_uniq_count.agg.txt", header = T, stringsAsFactors = F) %>% dplyr::select(KEGG_subcat, count)
group_II_III_common <- read.table("/media/chen/DATA/hasID_Cluster_II_III_common.agg.txt", header = T, stringsAsFactors = F) %>% dplyr::select(KEGG_subcat, count)

group_II_uniq.agg <- aggregate(count ~ KEGG_subcat, data = group_II_uniq, FUN = sum) %>% dplyr::mutate(freq = count/sum(group_II_uniq.agg$count))
group_III_uniq.agg <- aggregate(count ~ KEGG_subcat, data = group_III_uniq, FUN = sum) %>% dplyr::mutate(freq = count/sum(group_III_uniq.agg$count))
group_II_III_common.agg <- aggregate(count ~ KEGG_subcat, data = group_II_III_common, FUN = sum) %>% dplyr::mutate(freq = count/sum(group_II_III_common.agg$count))

tous_KEGG_BRITE_classification <- dplyr::bind_rows(group_II_uniq.agg, group_III_uniq.agg, group_II_III_common.agg) %>% dplyr::select(KEGG_subcat) %>% unique()

heatmap_matrix_machen <- function(df_2, df_3, df_2_3, df_classification) {
  df_2_p <- dplyr::inner_join(df_2, df_classification, by = "KEGG_subcat") %>% dplyr::select(KEGG_subcat, freq) %>%dplyr::rename(Group_II_uniq = freq)
  df_2_a <- dplyr::anti_join(df_classification, df_2, by = "KEGG_subcat") %>% dplyr::mutate(Group_II_uniq = 0)
  df_2_tous <- dplyr::bind_rows(df_2_p, df_2_a) %>% dplyr::arrange(desc(KEGG_subcat))
  
  df_3_p <- dplyr::inner_join(df_3, df_classification, by = "KEGG_subcat") %>% dplyr::select(KEGG_subcat, freq) %>%dplyr::rename(Group_III_uniq = freq)
  df_3_a <- dplyr::anti_join(df_classification, df_3, by = "KEGG_subcat") %>% dplyr::mutate(Group_III_uniq = 0)
  df_3_tous <- dplyr::bind_rows(df_3_p, df_3_a) %>% dplyr::arrange(desc(KEGG_subcat))
  
  df_2_3_p <- dplyr::inner_join(df_2_3, df_classification, by = "KEGG_subcat") %>% dplyr::select(KEGG_subcat, freq) %>%dplyr::rename(Group_II_III_common = freq)
  df_2_3_a <- dplyr::anti_join(df_classification, df_2_3, by = "KEGG_subcat") %>% dplyr::mutate(Group_II_III_common = 0) 
  df_2_3_tous <- dplyr::bind_rows(df_2_3_p, df_2_3_a) %>% dplyr::arrange(desc(KEGG_subcat))
  
  df <- dplyr::bind_cols(df_2_tous, df_3_tous, df_2_3_tous) %>% dplyr::select(KEGG_subcat...1, Group_II_uniq, Group_III_uniq, Group_II_III_common)
  row.names(df) <- df$KEGG_subcat...1
  df$KEGG_subcat...1 <- NULL
  df.mx <- data.matrix(df, rownames.force = T)
  
  return(df.mx)
}

df_figure_2e.mx <- heatmap_matrix_machen(group_II_uniq.agg, group_III_uniq.agg, group_II_III_common.agg, tous_KEGG_BRITE_classification)

Heatmap(df_figure_2e.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "", row_title = "KEGG BRITE classification", show_row_names = T, column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))

## version 0 ##
list_kegg_subcat_cluster_II_uniq <- hasID_Cluster_II_uniq_count.agg %>% dplyr::select(KEGG_subcat) %>% unique()
list_kegg_subcat_cluster_III_uniq <- hasID_Cluster_III_uniq_count.agg %>% dplyr::select(KEGG_subcat) %>% unique()
list_kegg_subcat_cluster_II_III_common <- hasID_Cluster_II_III_common.agg %>% dplyr::select(KEGG_subcat) %>% unique()


list_shared_kegg_subcat <- dplyr::inner_join(list_kegg_subcat_cluster_II_uniq, list_kegg_subcat_cluster_III_uniq, by = "KEGG_subcat")
list_shared_kegg_subcat <- dplyr::inner_join(list_shared_kegg_subcat, list_kegg_subcat_cluster_II_III_common, by = "KEGG_subcat")

hasID_Cluster_II_uniq_count.agg.sansID <- hasID_Cluster_II_uniq_count.agg %>% dplyr::select(KEGG_subcat, count)
hasID_Cluster_II_uniq_count.agg.sansID.agg <- aggregate(count ~ KEGG_subcat, data = hasID_Cluster_II_uniq_count.agg.sansID, FUN = sum)
hasID_Cluster_II_uniq_count.agg.sansID.agg.shared <- dplyr::inner_join(hasID_Cluster_II_uniq_count.agg.sansID.agg, list_shared_kegg_subcat, by = "KEGG_subcat") %>% rename(cluster_II_uniq = count) %>% arrange(desc(KEGG_subcat)) %>% mutate(cluster_II_uniq = cluster_II_uniq/dim(list_kegg_subcat_cluster_II_uniq)[1])

hasID_Cluster_III_uniq_count.agg.sansID <- hasID_Cluster_III_uniq_count.agg %>% dplyr::select(KEGG_subcat, count)
hasID_Cluster_III_uniq_count.agg.sansID.agg <- aggregate(count ~ KEGG_subcat, data = hasID_Cluster_III_uniq_count.agg.sansID, FUN = sum) 
hasID_Cluster_III_uniq_count.agg.sansID.agg.shared <- dplyr::inner_join(hasID_Cluster_III_uniq_count.agg.sansID.agg, list_shared_kegg_subcat, by = "KEGG_subcat") %>% rename(cluster_III_uniq = count) %>% arrange(desc(KEGG_subcat)) %>% mutate(cluster_III_uniq = cluster_III_uniq/dim(list_kegg_subcat_cluster_III_uniq)[1])

hasID_Cluster_II_III_common.agg.sansID <- hasID_Cluster_II_III_common.agg %>% dplyr::select(KEGG_subcat, count)
hasID_Cluster_II_III_common.agg.sansID.agg <- aggregate(count ~ KEGG_subcat, data = hasID_Cluster_II_III_common.agg.sansID, FUN = sum) %>% arrange(desc(KEGG_subcat))
hasID_Cluster_II_III_common.agg.sansID.agg.shared <- dplyr::inner_join(hasID_Cluster_II_III_common.agg.sansID.agg, list_shared_kegg_subcat, by = "KEGG_subcat") %>% rename(common = count) %>% arrange(desc(KEGG_subcat)) %>% mutate(common = common/dim(list_kegg_subcat_cluster_II_III_common)[1])

df_shared_pt_in_uniq_hsaID <- dplyr::bind_cols(hasID_Cluster_II_uniq_count.agg.sansID.agg.shared, hasID_Cluster_III_uniq_count.agg.sansID.agg.shared, hasID_Cluster_II_III_common.agg.sansID.agg.shared) %>% dplyr::select(KEGG_subcat...1, cluster_II_uniq, cluster_III_uniq, common)
row.names(df_shared_pt_in_uniq_hsaID) <- df_shared_pt_in_uniq_hsaID$KEGG_subcat...1
df_shared_pt_in_uniq_hsaID$KEGG_subcat...1 <- NULL
df_shared_pt_in_uniq_hsaID.mx <- data.matrix(df_shared_pt_in_uniq_hsaID, rownames.force = T)

write.table(df_shared_pt_in_uniq_hsaID.mx, file = "/media/chen/DATA/evoPath/df/C7/comparaison_cluster_2_3/df_shared_pt_in_uniq_hsaID.mx.txt", row.names = T, col.names = T, quote = F, sep = "\t")

Heatmap(df_shared_pt_in_uniq_hsaID.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.5), column_title = "clinical condition", row_title = "Cytokines", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
