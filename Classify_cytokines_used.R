cART_EC_cytokines_used <- read.table("/media/chen/DATA/evoPath/df/C7/c7_cART_EC_pt_cytokines.tsv", header = T, stringsAsFactors = F)
rownames(cART_EC_cytokines_used) <- cART_EC_cytokines_used$Cytokine
cART_EC_cytokines_used$Cytokine <- NULL
cART_EC_cytokines_used.mx <- data.matrix(cART_EC_cytokines_used, rownames.force = T)
  
# heatmap
pdf("/media/chen/DATA/evoPath/Abb/msigdb_cART_EC_pt_cytokines_used_v1.pdf", height = 4.5, width = 4)  
Heatmap(cART_EC_cytokines_used.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.5), column_title = "clinical condition", row_title = "Cytokines", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/msigdb_cART_EC_pt_cytokines_used_v1.svg", height = 4.5, width = 4)  
Heatmap(cART_EC_cytokines_used.mx, name = "Frequency", rect_gp = gpar(col = "black", lwd = 0.5), column_title = "clinical condition", row_title = "Cytokines", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
dev.off()
