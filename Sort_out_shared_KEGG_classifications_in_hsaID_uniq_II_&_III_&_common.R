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
