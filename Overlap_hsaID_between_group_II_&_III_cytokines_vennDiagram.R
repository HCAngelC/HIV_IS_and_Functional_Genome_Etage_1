#1. Plot venn diagram
tous_hasID_cluster_II <- dplyr::bind_rows(c7_kegg_short_CXCR5, c7_kegg_long_CXCR5, c7_kegg_short_IFNB, c7_kegg_long_IFNB, c7_kegg_short_IFNG, c7_kegg_long_IFNG, c7_kegg_short_IL10, c7_kegg_long_IL10, c7_kegg_short_TGFB, c7_kegg_long_TGFB) %>% dplyr::select(ID) %>% unique()

tous_hasID_cluster_III <- dplyr::bind_rows(c7_kegg_subcat_long_CXCL4, c7_kegg_subcat_long_IL18, c7_kegg_subcat_long_IL15, c7_kegg_subcat_long_IL7, c7_kegg_subcat_long_IL1, c7_kegg_subcat_long_IL2, c7_kegg_subcat_long_IFNA) %>% dplyr::select(ID) %>% unique()

overlap_hsaID_cluster2_3 <- dplyr::inner_join(tous_hasID_cluster_II, tous_hasID_cluster_III, by = "ID") %>% unique()

hasID_Cluster_II_uniq <- dplyr::anti_join(tous_hasID_cluster_II, overlap_hsaID_cluster2_3, by = "ID") %>% unique()
hasID_Cluster_III_uniq <- dplyr::anti_join(tous_hasID_cluster_III, overlap_hsaID_cluster2_3, by = "ID") %>% unique()

VennDiagram::draw.pairwise.venn(dim(tous_hasID_cluster_II)[1], dim(tous_hasID_cluster_III)[1], dim(overlap_hsaID_cluster2_3)[1], fill = c("#0073C2FF", "#CD534CFF"))
dev.off()
