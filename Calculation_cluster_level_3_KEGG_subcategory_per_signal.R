# 1. import df cluster 3 KEGG subcat
c7_kegg_subcat_long_CXCL4 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_CXCL4.csv", header = T, stringsAsFactors = F)
c7_kegg_subcat_long_IL18 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_IL18.csv", header = T, stringsAsFactors = F)
c7_kegg_subcat_long_IL15 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_IL15.csv", header = T, stringsAsFactors = F)
c7_kegg_subcat_long_IL7 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_IL7.csv", header = T, stringsAsFactors = F)
c7_kegg_subcat_long_IL1 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_IL1.csv", header = T, stringsAsFactors = F)
c7_kegg_subcat_long_IL2 <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_IL2.csv", header = T, stringsAsFactors = F)
c7_kegg_subcat_long_IFNA <- read.table("/media/chen/DATA/evoPath/df/C7/df_cluster_3/KEGG_dans_C7/kegg_c7_long_cluster_3_IFNA.csv", header = T, stringsAsFactors = F)

# 2. heatmap signal vs KEGG_subcat
##############################################################
# Author: Heng-Chang Chen
# Date: April 2022
##############################################################
# Input: dataframe containing each enriched hsa ID, KEGG subcategory & count number
# Object: generate the matrix to plot the clustering heatmap
##############################################################
# Custom R functions
**************************************************************

make_matrix_kegg_subcat_cluster_3 <- function(df, ref) {
  df <- df %>% dplyr::select(KEGG_subcat, count)
  df.agg <- aggregate(count ~ KEGG_subcat, data = df, FUN = sum)
  
  df.present <- dplyr::inner_join(df.agg, ref, by = "KEGG_subcat")
  df.absent <- dplyr::anti_join(ref, df.agg, by = "KEGG_subcat") %>% dplyr::mutate(count = 0)
  
  df.tous <- dplyr::bind_rows(df.present, df.absent) %>% arrange(desc(KEGG_subcat))
  
  return(df.tous)
}
**************************************************************
tous_term_KEGG_subcat_cluster_3 <- dplyr::bind_rows(c7_kegg_subcat_long_CXCL4, c7_kegg_subcat_long_IL18, c7_kegg_subcat_long_IL15, c7_kegg_subcat_long_IL7, c7_kegg_subcat_long_IL1, c7_kegg_subcat_long_IL2, c7_kegg_subcat_long_IFNA) %>% select(KEGG_subcat) %>% unique()

c7_kegg_tous_subcat_long_CXCL4 <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_CXCL4, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(CXCL4 = count)
c7_kegg_tous_subcat_long_IL18 <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_IL18, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(IL18 = count)
c7_kegg_tous_subcat_long_IL15 <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_IL15, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(IL15 = count)
c7_kegg_tous_subcat_long_IL7 <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_IL7, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(IL7 = count)
c7_kegg_tous_subcat_long_IL1 <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_IL1, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(IL1 = count)
c7_kegg_tous_subcat_long_IL2 <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_IL2, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(IL2 = count)
c7_kegg_tous_subcat_long_IFNA <- make_matrix_kegg_subcat_cluster_3(c7_kegg_subcat_long_IFNA, tous_term_KEGG_subcat_cluster_3) %>% dplyr::rename(IFNA = count)

tous_c7_kegg_subcat_cluster_3 <- dplyr::bind_cols(c7_kegg_tous_subcat_long_CXCL4, c7_kegg_tous_subcat_long_IL18, c7_kegg_tous_subcat_long_IL15, c7_kegg_tous_subcat_long_IL7, c7_kegg_tous_subcat_long_IL1, c7_kegg_tous_subcat_long_IL2, c7_kegg_tous_subcat_long_IFNA) %>% dplyr::select(KEGG_subcat...1, CXCL4, IL18, IL15, IL7, IL1, IL2, IFNA)

row.names(tous_c7_kegg_subcat_cluster_3) <- tous_c7_kegg_subcat_cluster_3$KEGG_subcat...1
tous_c7_kegg_subcat_cluster_3$KEGG_subcat...1 <- NULL
tous_c7_kegg_subcat_cluster_3.mx <- data.matrix(tous_c7_kegg_subcat_cluster_3, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/tous_c7_kegg_subcat_cluster_3.mx.pdf", height = 5, width = 6)  
Heatmap(tous_c7_kegg_subcat_cluster_3.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "signal", row_title = "KEGG subcategories", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/tous_c7_kegg_subcat_cluster_3.mx.svg", height = 5, width = 6)  
Heatmap(tous_c7_kegg_subcat_cluster_3.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "signal", row_title = "KEGG subcategories", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()
