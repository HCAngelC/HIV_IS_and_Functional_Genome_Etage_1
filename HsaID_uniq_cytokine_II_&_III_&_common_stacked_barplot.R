tous_hsaID_KEGG_subcat_cluster_2 <- dplyr::bind_rows(c7_kegg_short_CXCR5, c7_kegg_long_CXCR5, c7_kegg_short_IFNB, c7_kegg_long_IFNB, c7_kegg_short_IFNG, c7_kegg_long_IFNG, c7_kegg_short_IL10, c7_kegg_long_IL10, c7_kegg_short_TGFB, c7_kegg_long_TGFB) %>% dplyr::select(ID, KEGG_subcat) %>% dplyr::mutate(count = 1)

tous_hsaID_KEGG_subcat_cluster_3 <- dplyr::bind_rows(c7_kegg_subcat_long_CXCL4, c7_kegg_subcat_long_IL18, c7_kegg_subcat_long_IL15, c7_kegg_subcat_long_IL7, c7_kegg_subcat_long_IL1, c7_kegg_subcat_long_IL2, c7_kegg_subcat_long_IFNA) %>% dplyr::select(ID, KEGG_subcat) %>% dplyr::mutate(count = 1)

tous_hsaID_KEGG_subcat_cluster_2_3_pool <- dplyr::bind_rows(tous_hsaID_KEGG_subcat_cluster_2, tous_hsaID_KEGG_subcat_cluster_3)

# cluster II cytokines
hasID_Cluster_II_uniq_count <- dplyr::inner_join(tous_hsaID_KEGG_subcat_cluster_2, hasID_Cluster_II_uniq, by = "ID")
hasID_Cluster_II_uniq_count.agg <- aggregate(count ~ ID + KEGG_subcat, data = hasID_Cluster_II_uniq_count, FUN = sum)

ggplot(hasID_Cluster_II_uniq_count.agg, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(55))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# cluster III cytokines
hasID_Cluster_III_uniq_count <- dplyr::inner_join(tous_hsaID_KEGG_subcat_cluster_3, hasID_Cluster_III_uniq, by = "ID")
hasID_Cluster_III_uniq_count.agg <- aggregate(count ~ ID + KEGG_subcat, data = hasID_Cluster_III_uniq_count, FUN = sum)

ggplot(hasID_Cluster_III_uniq_count.agg, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Reds"))(35))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# shared hsaID
hasID_Cluster_II_III_common <- dplyr::inner_join(tous_hsaID_KEGG_subcat_cluster_2_3_pool, overlap_hsaID_cluster2_3, by = "ID")
hasID_Cluster_II_III_common.agg <- aggregate(count ~ ID + KEGG_subcat, data = hasID_Cluster_II_III_common, FUN = sum)

ggplot(hasID_Cluster_II_III_common.agg, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Greens"))(62))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
