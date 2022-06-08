#1. Import gene list in top 4 KEGG classification

list_gene_cancer_specific_type <- read.table("/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/Cancer_specific_type.csv", header = T, stringsAsFactors = T) %>% dplyr::select(Gene_ID, Gene_name) %>% unique()
list_gene_signal_transduction <- read.table("/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/Signal_transduction.csv", header = T, stringsAsFactors = T) %>% dplyr::select(Gene_ID, Gene_name) %>% unique()
list_gene_immune_system <- read.table("/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/Immune_system.csv", header = T, stringsAsFactors = T) %>% dplyr::select(Gene_ID, Gene_name) %>% unique()
list_gene_infectious_disease_viral <- read.table("/media/chen/DATA/evoPath/df/Top_4_KEGG_classifications/Infectious_disease_viral.csv", header = T, stringsAsFactors = T) %>% dplyr::select(Gene_ID, Gene_name) %>% unique()

#2. Plot venn diagram
overlap_1cancer_2singal <- merge(list_gene_cancer_specific_type, list_gene_signal_transduction, by = "Gene_name")
overlap_1cancer_3immune <- merge(list_gene_cancer_specific_type, list_gene_immune_system, by = "Gene_name")
overlap_2singal_3immune <- merge(list_gene_signal_transduction, list_gene_immune_system, by = "Gene_name")

n123 <- merge(overlap_1cancer_2singal, list_gene_immune_system, by = "Gene_name")

count_n123 <- dim(merge(overlap_1cancer_2singal, list_gene_immune_system, by = "Gene_name"))[1]
count_n124 <- dim(merge(overlap_1cancer_2singal, list_gene_infectious_disease_viral, by = "Gene_name"))[1]
count_n134 <- dim(merge(overlap_1cancer_3immune, list_gene_infectious_disease_viral, by = "Gene_name"))[1]
count_n234 <- dim(merge(overlap_2singal_3immune, list_gene_infectious_disease_viral, by = "Gene_name"))[1]

count_n12 <- dim(overlap_1cancer_2singal)[1]
count_n13 <- dim(overlap_1cancer_3immune)[1]
count_n14 <- dim(merge(list_gene_cancer_specific_type, list_gene_infectious_disease_viral, by = "Gene_name"))[1]
count_n23 <- dim(overlap_2singal_3immune)[1]
count_n24 <- dim(merge(list_gene_signal_transduction, list_gene_infectious_disease_viral, by = "Gene_name"))[1]
count_n34 <- dim(merge(list_gene_immune_system, list_gene_infectious_disease_viral, by = "Gene_name"))[1]

count_n1234 <- dim(merge(n123, list_gene_infectious_disease_viral, by = "Gene_name"))[1]

pdf("/media/chen/DATA/evoPath/Abb/Top_4_KEGG_pt_venn_diagram.pdf")  
draw.quad.venn(dim(list_gene_cancer_specific_type)[1], dim(list_gene_signal_transduction)[1], dim(list_gene_immune_system)[1], dim(list_gene_infectious_disease_viral)[1], count_n12, count_n13, count_n14, count_n23, count_n24, count_n34, count_n123, count_n124, count_n134, count_n234, count_n1234, category = c("Cancer specific type", "Signal transduction", "Immune system", "Infectious disease viral"), fill = c("#0073C2FF", "#CD534CFF", "#EFC000FF", "#868686FF"))
dev.off()
