# 1. Sort KEGG subcategory in Signal transduction, Immune system, Cancer specific types & Infection viral
calcu_pie_chart_signal_transduction <- function(df_c7_kegg_subcat, df_rel_genes) {
  df_bait <- df_c7_kegg_subcat %>% dplyr::filter(KEGG_subcat == "Signal_transduction") %>% dplyr::select(ID)
  row.names(df_rel_genes) <- NULL
  
  df_merge <- dplyr::inner_join(df_rel_genes, df_bait, by = "ID") %>% dplyr::select(Description) %>% dplyr::mutate(count = 1)
  
  df_merge.agg <- aggregate(count ~ Description, data = df_merge, FUN = sum)
  
  return(df_merge.agg)
}

calcu_pie_chart_Immune_system <- function(df_c7_kegg_subcat, df_rel_genes) {
  df_bait <- df_c7_kegg_subcat %>% dplyr::filter(KEGG_subcat == "Immune_system") %>% dplyr::select(ID)
  row.names(df_rel_genes) <- NULL
  
  df_merge <- dplyr::inner_join(df_rel_genes, df_bait, by = "ID") %>% dplyr::select(Description) %>% dplyr::mutate(count = 1)
  
  df_merge.agg <- aggregate(count ~ Description, data = df_merge, FUN = sum)
  
  return(df_merge.agg)
}

calcu_pie_chart_Cancer_specific_types <- function(df_c7_kegg_subcat, df_rel_genes) {
  df_bait <- df_c7_kegg_subcat %>% dplyr::filter(KEGG_subcat == "Cancer_specific_types") %>% dplyr::select(ID)
  row.names(df_rel_genes) <- NULL
  
  df_merge <- dplyr::inner_join(df_rel_genes, df_bait, by = "ID") %>% dplyr::select(Description) %>% dplyr::mutate(count = 1)
  
  df_merge.agg <- aggregate(count ~ Description, data = df_merge, FUN = sum)
  
  return(df_merge.agg)
}

calcu_pie_chart_Infectious_disease_viral <- function(df_c7_kegg_subcat, df_rel_genes) {
  df_bait <- df_c7_kegg_subcat %>% dplyr::filter(KEGG_subcat == "Infectious_disease_viral") %>% dplyr::select(ID)
  row.names(df_rel_genes) <- NULL
  
  df_merge <- dplyr::inner_join(df_rel_genes, df_bait, by = "ID") %>% dplyr::select(Description) %>% dplyr::mutate(count = 1)
  
  df_merge.agg <- aggregate(count ~ Description, data = df_merge, FUN = sum)
  
  return(df_merge.agg)
}

# Signal_transduction
Descrip_count_cART_short_CXCR5_st <- calcu_pie_chart_signal_transduction(c7_kegg_short_CXCR5, kegg_cART_short_CXCR5_rel_genes) %>% dplyr::mutate(period = "short", signal = "CXCR5")
Descrip_count_cART_long_CXCR5_st <- calcu_pie_chart_signal_transduction(c7_kegg_long_CXCR5, kegg_cART_long_CXCR5_rel_genes) %>% dplyr::mutate(period = "long", signal = "CXCR5")
#Descrip_count_cART_tous_CXCR5 <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5, Descrip_count_cART_long_CXCR5)

Descrip_count_cART_short_IFNB_st <- calcu_pie_chart_signal_transduction(c7_kegg_short_IFNB, kegg_cART_short_IFNB_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNB")
Descrip_count_cART_long_IFNB_st <- calcu_pie_chart_signal_transduction(c7_kegg_long_IFNB, kegg_cART_long_IFNB_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNB") 
#Descrip_count_cART_tous_IFNB <- dplyr::bind_rows(Descrip_count_cART_short_IFNB, Descrip_count_cART_long_IFNB) 

#Descrip_count_cART_short_IFNG_st <- calcu_pie_chart_signal_transduction(c7_kegg_short_IFNG, kegg_cART_short_IFNG_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNG") # not found
Descrip_count_cART_long_IFNG_st <- calcu_pie_chart_signal_transduction(c7_kegg_long_IFNG, kegg_cART_long_IFNG_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNG")
#Descrip_count_cART_tous_IFNG <- dplyr::bind_rows(Descrip_count_cART_short_IFNG, Descrip_count_cART_long_IFNG) 

Descrip_count_cART_short_IL10_st <- calcu_pie_chart_signal_transduction(c7_kegg_short_IL10, kegg_cART_short_IL10_rel_genes) %>% dplyr::mutate(period = "short", signal = "IL10")
Descrip_count_cART_long_IL10_st <- calcu_pie_chart_signal_transduction(c7_kegg_long_IL10, kegg_cART_long_IL10_rel_genes) %>% dplyr::mutate(period = "long", signal = "IL10")
#Descrip_count_cART_tous_IL10 <- dplyr::bind_rows(Descrip_count_cART_short_IL10, Descrip_count_cART_long_IL10) 

Descrip_count_cART_short_TGFB_st <- calcu_pie_chart_signal_transduction(c7_kegg_short_TGFB, kegg_cART_short_TGFB_rel_genes) %>% dplyr::mutate(period = "short", signal = "TGFB")
Descrip_count_cART_long_TGFB_st <- calcu_pie_chart_signal_transduction(c7_kegg_long_TGFB, kegg_cART_long_TGFB_rel_genes) %>% dplyr::mutate(period = "long", signal = "TGFB")
#Descrip_count_cART_tous_TGFB <- dplyr::bind_rows(Descrip_count_cART_short_TGFB, Descrip_count_cART_long_TGFB) 

# Immune_system
#Descrip_count_cART_short_CXCR5_is <- calcu_pie_chart_Immune_system(c7_kegg_short_CXCR5, kegg_cART_short_CXCR5_rel_genes) %>% dplyr::mutate(period = "short", signal = "CXCR5") # not found
Descrip_count_cART_long_CXCR5_is <- calcu_pie_chart_Immune_system(c7_kegg_long_CXCR5, kegg_cART_long_CXCR5_rel_genes) %>% dplyr::mutate(period = "long", signal = "CXCR5")
#Descrip_count_cART_tous_CXCR5 <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5, Descrip_count_cART_long_CXCR5)

Descrip_count_cART_short_IFNB_is <- calcu_pie_chart_Immune_system(c7_kegg_short_IFNB, kegg_cART_short_IFNB_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNB")
#Descrip_count_cART_long_IFNB_is <- calcu_pie_chart_Immune_system(c7_kegg_long_IFNB, kegg_cART_long_IFNB_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNB")  # not found
#Descrip_count_cART_tous_IFNB <- dplyr::bind_rows(Descrip_count_cART_short_IFNB, Descrip_count_cART_long_IFNB) 

Descrip_count_cART_short_IFNG_is <- calcu_pie_chart_Immune_system(c7_kegg_short_IFNG, kegg_cART_short_IFNG_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNG")
Descrip_count_cART_long_IFNG_is <- calcu_pie_chart_Immune_system(c7_kegg_long_IFNG, kegg_cART_long_IFNG_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNG")
#Descrip_count_cART_tous_IFNG <- dplyr::bind_rows(Descrip_count_cART_short_IFNG, Descrip_count_cART_long_IFNG) 

Descrip_count_cART_short_IL10_is <- calcu_pie_chart_Immune_system(c7_kegg_short_IL10, kegg_cART_short_IL10_rel_genes) %>% dplyr::mutate(period = "short", signal = "IL10")
Descrip_count_cART_long_IL10_is <- calcu_pie_chart_Immune_system(c7_kegg_long_IL10, kegg_cART_long_IL10_rel_genes) %>% dplyr::mutate(period = "long", signal = "IL10")
#Descrip_count_cART_tous_IL10 <- dplyr::bind_rows(Descrip_count_cART_short_IL10, Descrip_count_cART_long_IL10) 

Descrip_count_cART_short_TGFB_is <- calcu_pie_chart_Immune_system(c7_kegg_short_TGFB, kegg_cART_short_TGFB_rel_genes) %>% dplyr::mutate(period = "short", signal = "TGFB")
#Descrip_count_cART_long_TGFB_is <- calcu_pie_chart_Immune_system(c7_kegg_long_TGFB, kegg_cART_long_TGFB_rel_genes) %>% dplyr::mutate(period = "long", signal = "TGFB")  # not found
#Descrip_count_cART_tous_TGFB <- dplyr::bind_rows(Descrip_count_cART_short_TGFB, Descrip_count_cART_long_TGFB) 

# Cancer_specific_types
Descrip_count_cART_short_CXCR5_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_short_CXCR5, kegg_cART_short_CXCR5_rel_genes) %>% dplyr::mutate(period = "short", signal = "CXCR5")
#Descrip_count_cART_long_CXCR5_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_long_CXCR5, kegg_cART_long_CXCR5_rel_genes) %>% dplyr::mutate(period = "long", signal = "CXCR5")  # not found
#Descrip_count_cART_tous_CXCR5 <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5, Descrip_count_cART_long_CXCR5)

Descrip_count_cART_short_IFNB_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_short_IFNB, kegg_cART_short_IFNB_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNB")
Descrip_count_cART_long_IFNB_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_long_IFNB, kegg_cART_long_IFNB_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNB")
#Descrip_count_cART_tous_IFNB <- dplyr::bind_rows(Descrip_count_cART_short_IFNB, Descrip_count_cART_long_IFNB) 

#Descrip_count_cART_short_IFNG_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_short_IFNG, kegg_cART_short_IFNG_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNG") # not found
Descrip_count_cART_long_IFNG_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_long_IFNG, kegg_cART_long_IFNG_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNG")
#Descrip_count_cART_tous_IFNG <- dplyr::bind_rows(Descrip_count_cART_short_IFNG, Descrip_count_cART_long_IFNG) 

Descrip_count_cART_short_IL10_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_short_IL10, kegg_cART_short_IL10_rel_genes) %>% dplyr::mutate(period = "short", signal = "IL10")
#Descrip_count_cART_long_IL10_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_long_IL10, kegg_cART_long_IL10_rel_genes) %>% dplyr::mutate(period = "long", signal = "IL10") # not found
#Descrip_count_cART_tous_IL10 <- dplyr::bind_rows(Descrip_count_cART_short_IL10, Descrip_count_cART_long_IL10) 

#Descrip_count_cART_short_TGFB_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_short_TGFB, kegg_cART_short_TGFB_rel_genes) %>% dplyr::mutate(period = "short", signal = "TGFB") # not found
Descrip_count_cART_long_TGFB_ca <- calcu_pie_chart_Cancer_specific_types(c7_kegg_long_TGFB, kegg_cART_long_TGFB_rel_genes) %>% dplyr::mutate(period = "long", signal = "TGFB")
#Descrip_count_cART_tous_TGFB <- dplyr::bind_rows(Descrip_count_cART_short_TGFB, Descrip_count_cART_long_TGFB) 

# Infectious_disease_viral
Descrip_count_cART_short_CXCR5_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_short_CXCR5, kegg_cART_short_CXCR5_rel_genes) %>% dplyr::mutate(period = "short", signal = "CXCR5")
#Descrip_count_cART_long_CXCR5_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_long_CXCR5, kegg_cART_long_CXCR5_rel_genes) %>% dplyr::mutate(period = "long", signal = "CXCR5") # not found
#Descrip_count_cART_tous_CXCR5 <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5, Descrip_count_cART_long_CXCR5)

Descrip_count_cART_short_IFNB_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_short_IFNB, kegg_cART_short_IFNB_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNB")
Descrip_count_cART_long_IFNB_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_long_IFNB, kegg_cART_long_IFNB_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNB")
#Descrip_count_cART_tous_IFNB <- dplyr::bind_rows(Descrip_count_cART_short_IFNB, Descrip_count_cART_long_IFNB) 

#Descrip_count_cART_short_IFNG_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_short_IFNG, kegg_cART_short_IFNG_rel_genes) %>% dplyr::mutate(period = "short", signal = "IFNG") # not found
Descrip_count_cART_long_IFNG_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_long_IFNG, kegg_cART_long_IFNG_rel_genes) %>% dplyr::mutate(period = "long", signal = "IFNG")
#Descrip_count_cART_tous_IFNG <- dplyr::bind_rows(Descrip_count_cART_short_IFNG, Descrip_count_cART_long_IFNG) 

Descrip_count_cART_short_IL10_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_short_IL10, kegg_cART_short_IL10_rel_genes) %>% dplyr::mutate(period = "short", signal = "IL10")
Descrip_count_cART_long_IL10_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_long_IL10, kegg_cART_long_IL10_rel_genes) %>% dplyr::mutate(period = "long", signal = "IL10")
#Descrip_count_cART_tous_IL10 <- dplyr::bind_rows(Descrip_count_cART_short_IL10, Descrip_count_cART_long_IL10) 

#Descrip_count_cART_short_TGFB_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_short_TGFB, kegg_cART_short_TGFB_rel_genes) %>% dplyr::mutate(period = "short", signal = "TGFB") # not found
#Descrip_count_cART_long_TGFB_inf <- calcu_pie_chart_Infectious_disease_viral(c7_kegg_long_TGFB, kegg_cART_long_TGFB_rel_genes) %>% dplyr::mutate(period = "long", signal = "TGFB") # not found
#Descrip_count_cART_tous_TGFB <- dplyr::bind_rows(Descrip_count_cART_short_TGFB, Descrip_count_cART_long_TGFB) 

# 2. Plot heatmaps
# 2.a Combine the descriptions based on short- & long period
#Signal transduction
Descrip_sig_trans_short <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5_st, Descrip_count_cART_short_IFNB_st, Descrip_count_cART_short_IL10_st, Descrip_count_cART_short_TGFB_st) 
Descrip_sig_trans_long <- dplyr::bind_rows(Descrip_count_cART_long_CXCR5_st, Descrip_count_cART_long_IFNB_st, Descrip_count_cART_long_IFNG_st, Descrip_count_cART_long_IL10_st, Descrip_count_cART_long_TGFB_st)
Descript_sig_trans_items <- dplyr::bind_rows(Descrip_sig_trans_short, Descrip_sig_trans_long) %>% dplyr::select(Description) %>% unique()

#Immune system
Descrip_imm_short <- dplyr::bind_rows(Descrip_count_cART_short_IFNB_is, Descrip_count_cART_short_IFNG_is, Descrip_count_cART_short_IL10_is, Descrip_count_cART_short_TGFB_is) 
Descrip_imm_long <- dplyr::bind_rows(Descrip_count_cART_long_IFNG_is, Descrip_count_cART_long_IL10_is)
Descript_imm_items <- dplyr::bind_rows(Descrip_imm_short, Descrip_imm_long) %>% dplyr::select(Description) %>% unique()

#Cancer specific
Descrip_can_short <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5_ca, Descrip_count_cART_short_IFNB_ca, Descrip_count_cART_short_IL10_ca) 
Descrip_can_long <- dplyr::bind_rows(Descrip_count_cART_long_IFNB_ca, Descrip_count_cART_long_IFNG_ca, Descrip_count_cART_long_TGFB_ca)
Descript_can_items <- dplyr::bind_rows(Descrip_can_short, Descrip_can_long) %>% dplyr::select(Description) %>% unique()

#Infection viral
Descrip_inf_short <- dplyr::bind_rows(Descrip_count_cART_short_CXCR5_inf, Descrip_count_cART_short_IFNB_inf, Descrip_count_cART_short_IL10_inf) 
Descrip_inf_long <- dplyr::bind_rows(Descrip_count_cART_long_IFNB_inf, Descrip_count_cART_long_IFNG_inf, Descrip_count_cART_long_IL10_inf)
Descript_inf_items <- dplyr::bind_rows(Descrip_inf_short, Descrip_inf_long) %>% dplyr::select(Description) %>% unique()

# 2.b Function for preparing the matrix
Description_pour_heatmap <- function(df_des_short, df_des_long, df_term) {
  df_des_short <- df_des_short %>% dplyr::select(Description, count)
  df_des_short.agg <- aggregate(count ~ Description, data = df_des_short, FUN = sum) %>% dplyr::rename(short = count)
  df_des_short.sans <- dplyr::anti_join(df_term, df_des_short.agg) %>% dplyr::mutate(short = 0)
  df_des_short.tous <- dplyr::bind_rows(df_des_short.agg , df_des_short.sans) %>% dplyr::arrange(desc(Description))
  
  df_des_long <- df_des_long %>% dplyr::select(Description, count)
  df_des_long.agg <- aggregate(count ~ Description, data = df_des_long, FUN = sum) %>% dplyr::rename(long = count)
  df_des_long.sans <- dplyr::anti_join(df_term, df_des_long.agg) %>% dplyr::mutate(long = 0)
  df_des_long.tous <- dplyr::bind_rows(df_des_long.agg , df_des_long.sans) %>% dplyr::arrange(desc(Description))
  
  df_tous <- dplyr::bind_cols(df_des_short.tous, df_des_long.tous)
  df_tous <- df_tous %>% dplyr::select(Description...1, short, long)
  row.names(df_tous) <- df_tous$Description...1
  df_tous$Description...1 <- NULL
  df_tous.mx <- data.matrix(df_tous, rownames.force = T)
  
  return(df_tous.mx)
}

Heatmap_description_signal_transduction <- Description_pour_heatmap(Descrip_sig_trans_short, Descrip_sig_trans_long, Descript_sig_trans_items)
pdf("/media/chen/DATA/evoPath/Abb/Heatmap_description_signal_transduction.pdf", height = 3, width = 5.4)  
Heatmap(Heatmap_description_signal_transduction, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Signal transduction", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()
svg("/media/chen/DATA/evoPath/Abb/Heatmap_description_signal_transduction.svg", height = 3, width = 5.4)  
Heatmap(Heatmap_description_signal_transduction, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Signal transduction", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()

Heatmap_description_immune_system <- Description_pour_heatmap(Descrip_imm_short, Descrip_imm_long, Descript_imm_items)
pdf("/media/chen/DATA/evoPath/Abb/Heatmap_description_immune_system.pdf", height = 2.5, width = 5.5)  
Heatmap(Heatmap_description_immune_system, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Immune system", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()
svg("/media/chen/DATA/evoPath/Abb/Heatmap_description_immune_system.svg", height = 2.5, width = 5.5)  
Heatmap(Heatmap_description_immune_system, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Immune system", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()

Heatmap_description_cancer <- Description_pour_heatmap(Descrip_can_short, Descrip_can_long, Descript_can_items)
pdf("/media/chen/DATA/evoPath/Abb/Heatmap_description_cancer.pdf", height = 4, width = 5)  
Heatmap(Heatmap_description_cancer, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Cancer specific types", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()
svg("/media/chen/DATA/evoPath/Abb/Heatmap_description_cancer.svg", height = 4, width = 5)  
Heatmap(Heatmap_description_cancer, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Cancer specific types", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()

Heatmap_description_infection_viral <- Description_pour_heatmap(Descrip_inf_short, Descrip_inf_long, Descript_inf_items)
pdf("/media/chen/DATA/evoPath/Abb/Heatmap_description_infection_viral.pdf", height = 2.5, width = 5.5)  
Heatmap(Heatmap_description_infection_viral, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Viral Infection", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()
svg("/media/chen/DATA/evoPath/Abb/Heatmap_description_infection_viral.svg", height = 2.5, width = 5.5)  
Heatmap(Heatmap_description_infection_viral, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.1), column_title = "Period", row_title = "Viral infection", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), show_column_dend = FALSE, column_names_gp = gpar(fontsize = 10))
dev.off()
