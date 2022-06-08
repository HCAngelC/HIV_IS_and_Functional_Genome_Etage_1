##############################################################
# Author: Heng-Chang Chen
# Date: June 2022
##############################################################
# Input: dataframes: lsit of genes in each top 4 KEGG classification, list of HIV-1-targeted genes in cART (untreated), cART (short), cART (long) and long-term EC.
# Object: Heatmap representation of the genes involved in top 4 KEGG classifications in each patient type
##############################################################
# Custom R functions
**************************************************************
fish_gens_dans_les_patients <- function(df_input, df_cART_ut, df_cART_st, df_cART_lt, df_ec) {
  df_input <- df_input %>% dplyr::select(Gene_name) %>% dplyr::rename(V1 = Gene_name) %>% unique()
  
  ut_p <- dplyr::inner_join(df_cART_ut, df_input, by = "V1") %>% dplyr::mutate(ut = "2")
  ut_a <- dplyr::anti_join(df_input, df_cART_ut, by = "V1") %>% dplyr::mutate(ut = "1")
  ut_tous <- dplyr::bind_rows(ut_p, ut_a) %>% dplyr::arrange(desc(V1))
  
  st_p <- dplyr::inner_join(df_cART_st, df_input, by = "V1") %>% dplyr::mutate(st = "2")
  st_a <- dplyr::anti_join(df_input, df_cART_st, by = "V1") %>% dplyr::mutate(st = "1")
  st_tous <- dplyr::bind_rows(st_p, st_a) %>% dplyr::arrange(desc(V1))
  
  lt_p <- dplyr::inner_join(df_cART_lt, df_input, by = "V1") %>% dplyr::mutate(lt = "2")
  lt_a <- dplyr::anti_join(df_input, df_cART_lt, by = "V1") %>% dplyr::mutate(lt = "1")
  lt_tous <- dplyr::bind_rows(lt_p, lt_a) %>% dplyr::arrange(desc(V1))
  
  names(df_ec) <- c("V1")
  
  ec_p <- dplyr::inner_join(df_ec, df_input, by = "V1") %>% dplyr::mutate(ec = "2")
  ec_a <- dplyr::anti_join(df_input, df_ec, by = "V1") %>% dplyr::mutate(ec = "1")
  ec_tous <- dplyr::bind_rows(ec_p, ec_a) %>% dplyr::arrange(desc(V1))
  
  df_tous <- dplyr::bind_cols(ut_tous, st_tous, lt_tous, ec_tous) %>% dplyr::select(V1...1, ut, st, lt, ec)

  row.names(df_tous) <- df_tous$V1...1
  df_tous$V1...1 <- NULL
  df_tous.mx <- data.matrix(df_tous, rownames.force = T)
  
  return(df_tous.mx)
}

**************************************************************

gene_present_cancer_specific_type <- fish_gens_dans_les_patients(list_gene_cancer_specific_type, cART_untreat, cART_short, cART_long, long_term_EC)
gene_present_signal_transduction <- fish_gens_dans_les_patients(list_gene_signal_transduction, cART_untreat, cART_short, cART_long, long_term_EC)
gene_present_immune_system <- fish_gens_dans_les_patients(list_gene_immune_system, cART_untreat, cART_short, cART_long, long_term_EC)
gene_present_infectious_disease_viral <- fish_gens_dans_les_patients(list_gene_infectious_disease_viral, cART_untreat, cART_short, cART_long, long_term_EC)

pdf("/media/chen/DATA/evoPath/Abb/gene_present_cancer_specific_type.mx.pdf", height = 3, width = 5)  
colors = structure(1:2, names = c("1", "2"))
Heatmap(gene_present_cancer_specific_type, name = "Cancer specific type", col = colors, rect_gp = gpar(col = "white", lwd = 0.5), column_title = "", row_title = "Genes", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

pdf("/media/chen/DATA/evoPath/Abb/gene_present_signal_transduction.mx.pdf", height = 3, width = 4.8)  
colors = structure(1:2, names = c("1", "2"))
Heatmap(gene_present_signal_transduction, name = "Signal_transduction", col = colors, rect_gp = gpar(col = "white", lwd = 0.5), column_title = "", row_title = "Genes", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

pdf("/media/chen/DATA/evoPath/Abb/gene_present_immune_system.mx.pdf", height = 3, width = 4.8)  
colors = structure(1:2, names = c("1", "2"))
Heatmap(gene_present_immune_system, name = "Immune_system", col = colors, rect_gp = gpar(col = "white", lwd = 0.5), column_title = "", row_title = "Genes", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

pdf("/media/chen/DATA/evoPath/Abb/gene_present_infectious_disease_viral.mx.pdf", height = 3, width = 5.2)  
colors = structure(1:2, names = c("1", "2"))
Heatmap(gene_present_immune_system, name = "Infectious_disease_viral", col = colors, rect_gp = gpar(col = "white", lwd = 0.5), column_title = "", row_title = "Genes", column_title_side = "bottom", row_dend_width  = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()
