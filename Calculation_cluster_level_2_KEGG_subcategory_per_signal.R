##############################################################
# Author: Heng-Chang Chen
# Date: March 2022
##############################################################
# Input: 
# Object:  
# ->
##############################################################
# Custom R functions
**************************************************************
calculation_KEGG_subcategory_pour_matrix <- function(df, ref) {
  df <- df %>% dplyr::select(KEGG_subcat, count)
  df.agg <- aggregate(count ~ KEGG_subcat, data = df, FUN = sum)
  df.agg.avoir <- dplyr::inner_join(df.agg, ref, by = "KEGG_subcat")
  df.sans <- dplyr::anti_join(ref, df.agg, by ="KEGG_subcat") %>% dplyr::mutate(count = 0)
  df.tous <- dplyr::bind_rows(df.agg.avoir, df.sans) %>% dplyr::arrange(desc(KEGG_subcat))
  
  return(df.tous)
}
****************************************************************
tous_term_KEGG_subcat <- dplyr::bind_rows(c7_kegg_short_TGFB, c7_kegg_long_TGFB, c7_kegg_short_IFNB, c7_kegg_long_IFNB, c7_kegg_short_IFNG, c7_kegg_long_IFNG, c7_kegg_short_IL10, c7_kegg_long_IL10, c7_kegg_short_IFNB, c7_kegg_long_FOXP3) %>% dplyr::select(KEGG_subcat) %>% unique()


calcu_KEGG_subcat_TGFB_short <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_short_TGFB, tous_term_KEGG_subcat) %>% dplyr::rename(TGFB_short = count)
calcu_KEGG_subcat_TGFB_long <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_long_TGFB, tous_term_KEGG_subcat) %>% dplyr::rename(TGFB_long = count)

calcu_KEGG_subcat_IFNB_short <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_short_IFNB, tous_term_KEGG_subcat) %>% dplyr::rename(IFNB_short = count)
calcu_KEGG_subcat_IFNB_long <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_long_IFNB, tous_term_KEGG_subcat) %>% dplyr::rename(IFNB_long = count)

calcu_KEGG_subcat_IFNG_short <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_short_IFNG, tous_term_KEGG_subcat) %>% dplyr::rename(IFNG_short = count)
calcu_KEGG_subcat_IFNG_long <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_long_IFNG, tous_term_KEGG_subcat) %>% dplyr::rename(IFNG_long = count)

calcu_KEGG_subcat_IL10_short <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_short_IL10, tous_term_KEGG_subcat) %>% dplyr::rename(IL10_short = count)
calcu_KEGG_subcat_IL10_long <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_long_IL10, tous_term_KEGG_subcat) %>% dplyr::rename(IL10_long = count)

calcu_KEGG_subcat_CXCR5_short <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_short_CXCR5, tous_term_KEGG_subcat) %>% dplyr::rename(CXCR5_short = count)
calcu_KEGG_subcat_CXCR5_long <- calculation_KEGG_subcategory_pour_matrix(c7_kegg_long_CXCR5, tous_term_KEGG_subcat) %>% dplyr::rename(CXCR5_long = count)

df_calcu_KEGG_subcat_tous_short_long <- dplyr::bind_cols(calcu_KEGG_subcat_TGFB_short, calcu_KEGG_subcat_TGFB_long, calcu_KEGG_subcat_IFNB_short, calcu_KEGG_subcat_IFNB_long, calcu_KEGG_subcat_IFNG_short, calcu_KEGG_subcat_IFNG_long, calcu_KEGG_subcat_IL10_short, calcu_KEGG_subcat_IL10_long, calcu_KEGG_subcat_CXCR5_short, calcu_KEGG_subcat_CXCR5_long) %>% dplyr::select(KEGG_subcat...1, TGFB_short, TGFB_long, IFNB_short, IFNB_long, IFNG_short, IFNG_long, IL10_short, IL10_long, CXCR5_short, CXCR5_long)

row.names(df_calcu_KEGG_subcat_tous_short_long) <- df_calcu_KEGG_subcat_tous_short_long$KEGG_subcat...1
df_calcu_KEGG_subcat_tous_short_long$KEGG_subcat...1 <- NULL

df_calcu_KEGG_subcat_tous_short_long.mx <- data.matrix(df_calcu_KEGG_subcat_tous_short_long, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/df_calcu_KEGG_subcat_tous_short_long.mx.pdf", height = 5.5, width = 7)  
#col_c7_kegg_subcat <- colorRamp2(c(0, 10, 20, 30, 40), c("Grey", "Blue", "Green", "#CCCC99","Yellow"))
Heatmap(df_calcu_KEGG_subcat_tous_short_long.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "signal", row_title = "KEGG subcategories", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()
