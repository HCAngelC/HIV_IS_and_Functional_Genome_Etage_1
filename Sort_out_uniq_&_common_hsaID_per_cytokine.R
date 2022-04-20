#1. compare the proportion of KEGG subcat from uniq. hsa ID terms per cytokine
# Common hsaID
tous_hsaID_uniq <- dplyr::bind_rows(c7_kegg_subcat_long_CXCL4.uniq, c7_kegg_subcat_long_IL1.uniq, c7_kegg_subcat_long_IL2.uniq, c7_kegg_subcat_long_IL7.uniq, c7_kegg_subcat_long_IL18.uniq, c7_kegg_subcat_long_IL15.uniq) %>% dplyr::select(ID, KEGG_subcat) %>% unique()

tous_hsaID <- dplyr::bind_rows(c7_kegg_subcat_long_CXCL4, c7_kegg_subcat_long_IL18, c7_kegg_subcat_long_IL15, c7_kegg_subcat_long_IL7, c7_kegg_subcat_long_IL1, c7_kegg_subcat_long_IL2, c7_kegg_subcat_long_IFNA) %>% dplyr::select(ID, KEGG_subcat) %>% unique()

tous_hsaID_rest <- dplyr::anti_join(tous_hsaID, tous_hsaID_uniq, by = "ID") %>% dplyr::select(KEGG_subcat) %>% dplyr::mutate(count = 1)
tous_hsaID_rest.agg <- aggregate(count ~ KEGG_subcat, data = tous_hsaID_rest, FUN = sum)

# Unique hsaID
c7_kegg_subcat_long_CXCL4.uniq <- c7_kegg_subcat_long_CXCL4 %>% dplyr::filter(ID == "hsa00562" | ID == "hsa03450" | ID == "hsa04012" | ID == "hsa04062" | ID == "hsa04070" | ID == "hsa04917" | ID == "hsa04933" | ID == "hsa04935" | ID == "hsa05207" | ID == "hsa05220" | ID == "hsa05223" | ID == "hsa05221") 
c7_kegg_subcat_long_CXCL4.uniq.agg <- aggregate(count ~ KEGG_subcat, data = c7_kegg_subcat_long_CXCL4.uniq, FUN = sum)

c7_kegg_subcat_long_IL1.uniq <- c7_kegg_subcat_long_IL1 %>% dplyr::filter(ID == "hsa03015" | ID == "hsa04340" | ID == "hsa04510" | ID == "hsa04512" | ID == "hsa04514" | ID == "hsa04640" | ID == "hsa05134" | ID == "hsa05222" | ID == "hsa05410" | ID == "hsa05414" | ID == "hsa05412") %>% dplyr::select(KEGG_subcat, count) 
c7_kegg_subcat_long_IL1.uniq.agg <- aggregate(count ~ KEGG_subcat, data = c7_kegg_subcat_long_IL1.uniq, FUN = sum)

c7_kegg_subcat_long_IL2.uniq <- c7_kegg_subcat_long_IL2 %>% dplyr::filter(ID == "hsa02010" | ID == "hsa03013" | ID == "hsa04015" | ID == "hsa04145" | ID == "hsa04390" | ID == "hsa04612" | ID == "hsa04530") %>% dplyr::select(KEGG_subcat, count) 
c7_kegg_subcat_long_IL2.uniq.agg <- aggregate(count ~ KEGG_subcat, data = c7_kegg_subcat_long_IL2.uniq, FUN = sum)

c7_kegg_subcat_long_IL7.uniq <- c7_kegg_subcat_long_IL7 %>% dplyr::filter(ID == "hsa03022" | ID == "hsa00604") %>% dplyr::select(KEGG_subcat, count)
c7_kegg_subcat_long_IL7.uniq.agg <- aggregate(count ~ KEGG_subcat, data = c7_kegg_subcat_long_IL7.uniq, FUN = sum)

c7_kegg_subcat_long_IL18.uniq <- c7_kegg_subcat_long_IL18 %>% dplyr::filter(ID == "hsa03040") %>% dplyr::select(KEGG_subcat, count)
 
c7_kegg_subcat_long_IL15.uniq <- c7_kegg_subcat_long_IL15 %>% dplyr::filter(ID == "hsa04072" | ID == "hsa04144" | ID == "hsa05131" | ID == "hsa05168" | ID == "hsa05132")  %>% dplyr::select(KEGG_subcat, count)
c7_kegg_subcat_long_IL15.uniq.agg <- aggregate(count ~ KEGG_subcat, data = c7_kegg_subcat_long_IL15.uniq, FUN = sum)

# Make a matrix
tous_kegg_subcat_uniq_common <- dplyr::bind_rows(c7_kegg_subcat_long_CXCL4, c7_kegg_subcat_long_IL18, c7_kegg_subcat_long_IL15, c7_kegg_subcat_long_IL7, c7_kegg_subcat_long_IL1, c7_kegg_subcat_long_IL2, c7_kegg_subcat_long_IFNA) %>% dplyr::select(KEGG_subcat) %>% unique()

make_matrix_kegg_subcat_uniq_common <- function(df, ref) {
  #df.present <- dplyr::inner_join(df, ref, by = "KEGG_subcat")
  df.absent <- dplyr::anti_join(ref, df, by = "KEGG_subcat") %>% dplyr::mutate(count = 0)
  
  df.tous <- dplyr::bind_rows(df, df.absent) %>% arrange(desc(KEGG_subcat))
  
  return(df.tous)
}

c7_kegg_subcat_CXCL4.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(c7_kegg_subcat_long_CXCL4.uniq.agg, tous_kegg_subcat_uniq_common) %>% dplyr::rename(CXCL4_uniq = count)
c7_kegg_subcat_IL1.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(c7_kegg_subcat_long_IL1.uniq.agg, tous_kegg_subcat_uniq_common) %>% dplyr::rename(IL1_uniq = count)
c7_kegg_subcat_IL2.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(c7_kegg_subcat_long_IL2.uniq.agg, tous_kegg_subcat_uniq_common) %>% dplyr::rename(IL2_uniq = count)
c7_kegg_subcat_IL7.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(c7_kegg_subcat_long_IL7.uniq.agg, tous_kegg_subcat_uniq_common) %>% dplyr::rename(IL7_uniq = count)
c7_kegg_subcat_IL15.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(c7_kegg_subcat_long_IL15.uniq.agg, tous_kegg_subcat_uniq_common) %>% dplyr::rename(IL15_uniq = count)
c7_kegg_subcat_IL18.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(c7_kegg_subcat_long_IL18.uniq, tous_kegg_subcat_uniq_common) %>% dplyr::rename(IL18_uniq = count)
c7_kegg_subcat_common.uniq_dans_tous <- make_matrix_kegg_subcat_uniq_common(tous_hsaID_rest.agg , tous_kegg_subcat_uniq_common) %>% dplyr::rename(Common = count)

df_tous_kegg_subcat_cluster_3_uniq_common <- dplyr::bind_cols(c7_kegg_subcat_CXCL4.uniq_dans_tous, c7_kegg_subcat_IL1.uniq_dans_tous, c7_kegg_subcat_IL2.uniq_dans_tous, c7_kegg_subcat_IL7.uniq_dans_tous, c7_kegg_subcat_IL15.uniq_dans_tous, c7_kegg_subcat_IL18.uniq_dans_tous, c7_kegg_subcat_common.uniq_dans_tous) %>% dplyr::select(KEGG_subcat...1, CXCL4_uniq, IL1_uniq, IL2_uniq, IL7_uniq, IL15_uniq, IL18_uniq, Common)
row.names(df_tous_kegg_subcat_cluster_3_uniq_common) <- df_tous_kegg_subcat_cluster_3_uniq_common$KEGG_subcat...1
df_tous_kegg_subcat_cluster_3_uniq_common$KEGG_subcat...1 <- NULL
df_tous_kegg_subcat_cluster_3_uniq_common.mx <- data.matrix(df_tous_kegg_subcat_cluster_3_uniq_common, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/df_tous_kegg_subcat_cluster_3_uniq_common.mx.pdf", height = 5, width = 6)  
Heatmap(df_tous_kegg_subcat_cluster_3_uniq_common.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "signal", row_title = "KEGG subcategories", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/df_tous_kegg_subcat_cluster_3_uniq_common.mx.svg", height = 5, width = 6)  
Heatmap(df_tous_kegg_subcat_cluster_3_uniq_common.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "signal", row_title = "KEGG subcategories", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()
