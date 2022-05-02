# 1. C7 immunological terms enrichment
m_df <- msigdbr(species = "Homo sapiens")

m_tC7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% dplyr::select(gs_name, entrez_gene)

msigdb_c7_cART_untreat <- data.frame(enricher(cART_untreat.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_cART_short <- data.frame(enricher(cART_short.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_cART_long <- data.frame(enricher(cART_long.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_EC_long <- data.frame(enricher(long_term_EC.enID$ENTREZID, TERM2GENE = m_tC7))

msigdb_c7_cART_untreat.key <- msigdb_c7_cART_untreat %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "untreat")
msigdb_c7_cART_untreat.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_untreat.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_untreat.key$GeneRatio, perl=T))
msigdb_c7_cART_untreat.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_untreat.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_untreat.key$BgRatio, perl=T))

msigdb_c7_cART_short.key <- msigdb_c7_cART_short %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "short")
msigdb_c7_cART_short.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$GeneRatio, perl=T))
msigdb_c7_cART_short.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$BgRatio, perl=T))

msigdb_c7_cART_long.key <- msigdb_c7_cART_long %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "long")
msigdb_c7_cART_long.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$GeneRatio, perl=T))
msigdb_c7_cART_long.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$BgRatio, perl=T))

msigdb_c7_EC_long.key <- msigdb_c7_EC_long %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "EC long")
msigdb_c7_EC_long.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_EC_long.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_EC_long.key$GeneRatio, perl=T))
msigdb_c7_EC_long.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_EC_long.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_EC_long.key$BgRatio, perl=T))

msigdb_cART_EC_tous <- dplyr::bind_rows(msigdb_c7_cART_untreat.key, msigdb_c7_cART_long.key, msigdb_c7_cART_short.key, msigdb_c7_EC_long.key) %>% dplyr::mutate(Ratio = GeneRatio/BgRatio)
msigdb_cART_EC_tous.cART_untreat <- dplyr::filter(msigdb_cART_EC_tous, cond == "cART untreat")
msigdb_cART_EC_tous.cART_long <- dplyr::filter(msigdb_cART_EC_tous, cond == "cART long")
msigdb_cART_EC_tous.cART_short <- dplyr::filter(msigdb_cART_EC_tous, cond == "cART short")
msigdb_cART_EC_tous.EC_long <- dplyr::filter(msigdb_cART_EC_tous, cond == "EC long")

msigdb_cART_EC_tous_ID_unique <- msigdb_cART_EC_tous %>% dplyr::select(ID) %>% unique()

msigdb_cART_EC_tous.cART_untreat_avoir <- dplyr::inner_join(msigdb_cART_EC_tous.cART_untreat, msigdb_cART_EC_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_EC_tous.cART_untreat_sans <- dplyr::anti_join(msigdb_cART_EC_tous_ID_unique, msigdb_cART_EC_tous.cART_untreat, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_EC_tous.cART_untreat.ID.Ratio <- dplyr::bind_rows(msigdb_cART_EC_tous.cART_untreat_avoir, msigdb_cART_EC_tous.cART_untreat_sans) %>% dplyr::rename(cART_untreat = Ratio) %>% arrange(desc(ID))
rownames(msigdb_cART_EC_tous.cART_untreat.ID.Ratio) <- NULL

msigdb_cART_EC_tous.cART_long_avoir <- dplyr::inner_join(msigdb_cART_EC_tous.cART_long, msigdb_cART_EC_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_EC_tous.cART_long_sans <- dplyr::anti_join(msigdb_cART_EC_tous_ID_unique, msigdb_cART_EC_tous.cART_long, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_EC_tous.cART_long.ID.Ratio <- dplyr::bind_rows(msigdb_cART_EC_tous.cART_long_avoir, msigdb_cART_EC_tous.cART_long_sans) %>% dplyr::rename(cART_long = Ratio) %>% arrange(desc(ID))
rownames(msigdb_cART_EC_tous.cART_long.ID.Ratio) <- NULL

msigdb_cART_EC_tous.cART_short_avoir <- dplyr::inner_join(msigdb_cART_EC_tous.cART_short, msigdb_cART_EC_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_EC_tous.cART_short_sans <- dplyr::anti_join(msigdb_cART_EC_tous_ID_unique, msigdb_cART_EC_tous.cART_short, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_EC_tous.cART_short.ID.Ratio <- dplyr::bind_rows(msigdb_cART_EC_tous.cART_short_avoir, msigdb_cART_EC_tous.cART_short_sans) %>% dplyr::rename(cART_short = Ratio) %>% arrange(desc(ID))
rownames(msigdb_cART_EC_tous.cART_short.ID.Ratio) <- NULL

msigdb_cART_EC_tous.EC_long_avoir <- dplyr::inner_join(msigdb_cART_EC_tous.EC_long, msigdb_cART_EC_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_EC_tous.EC_long_sans <- dplyr::anti_join(msigdb_cART_EC_tous_ID_unique, msigdb_cART_EC_tous.EC_long, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_EC_tous.EC_long.ID.Ratio <- dplyr::bind_rows(msigdb_cART_EC_tous.EC_long_avoir, msigdb_cART_EC_tous.EC_long_sans) %>% dplyr::rename(EC_long = Ratio) %>% arrange(desc(ID))
rownames(msigdb_cART_EC_tous.EC_long.ID.Ratio) <- NULL

msigdb_cART_EC_tous.ID.Ratio <- dplyr::bind_cols(msigdb_cART_EC_tous.cART_untreat.ID.Ratio, msigdb_cART_EC_tous.cART_short.ID.Ratio, msigdb_cART_EC_tous.cART_long.ID.Ratio, msigdb_cART_EC_tous.EC_long.ID.Ratio) %>% dplyr::select(ID...1, cART_untreat, cART_short, cART_long, EC_long)
row.names(msigdb_cART_EC_tous.ID.Ratio) <- msigdb_cART_EC_tous.ID.Ratio$ID...1
msigdb_cART_EC_tous.ID.Ratio$ID...1 <- NULL
msigdb_cART_EC_tous.ID.Ratio.mx <- data.matrix(msigdb_cART_EC_tous.ID.Ratio, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/msigdb_cART_EC_evolution_pathway.pdf", height = 5, width = 4)  
Heatmap(msigdb_cART_EC_tous.ID.Ratio.mx, name = "GeneRatio/BgRatio", rect_gp = gpar(col = "black", lwd = 0.01), column_title = "clinical condition", row_title = "C7: immunologic signatures", column_title_side = "bottom", row_dend_width = unit(4, "cm"), show_row_names = F, column_names_gp = gpar(fontsize = 10))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/msigdb_cART_EC_evolution_pathway.svg", height = 5, width = 4)  
Heatmap(msigdb_cART_EC_tous.ID.Ratio.mx, name = "GeneRatio/BgRatio", rect_gp = gpar(col = "black", lwd = 0.01), column_title = "clinical condition", row_title = "C7: immunologic signatures", column_title_side = "bottom", row_dend_width = unit(4, "cm"), show_row_names = F, column_names_gp = gpar(fontsize = 10))
dev.off()

# 2. Pathway overlapped
overlap_cART_c7_short_long <- merge(msigdb_c7_cART_short.key, msigdb_c7_cART_long.key, by = "ID")
count_n123 <- dim(merge(overlap_cART_c7_short_long, msigdb_c7_cART_untreat.key, by = "ID"))[1]

count_n12 <- dim(merge(msigdb_c7_cART_untreat.key, msigdb_c7_cART_short.key, by = "ID"))[1]
count_n23 <- dim(merge(msigdb_c7_cART_short.key, msigdb_c7_cART_long.key, by = "ID"))[1]
count_n13 <- dim(merge(msigdb_c7_cART_untreat.key, msigdb_c7_cART_long.key, by = "ID"))[1]

pdf("/media/chen/DATA/evoPath/Abb/msigdb_cART_c7_venn_diagram.pdf")  
draw.triple.venn(dim(msigdb_c7_cART_untreat.key)[1], dim(msigdb_c7_cART_short.key)[1], dim(msigdb_c7_cART_long.key)[1], count_n12, count_n23, count_n13, count_n123, category = c("untreat", "short", "long"), rotation = 1, reverse = FALSE, euler.d =
TRUE, scaled = TRUE, lwd = rep(2, 3), lty =
rep("solid", 3), col = rep("black", 3), fill = NULL,
alpha = rep(0.5, 3), label.col = rep("black", 7))
dev.off()

***include EC***
list_ID <- list(cART_untreat = msigdb_c7_cART_untreat.key$ID,
                cART_short = msigdb_c7_cART_short.key$ID,
                cART_long = msigdb_c7_cART_long.key$ID,
                EC_long = msigdb_c7_EC_long.key$ID)

overlap_cART_c7_short_long <- merge(msigdb_c7_cART_short.key, msigdb_c7_cART_long.key, by = "ID")
overlap_cART_c7_untreat_short <- merge(msigdb_c7_cART_untreat.key, msigdb_c7_cART_short.key, by = "ID")
overlap_cART_c7_untreat_long <- merge(msigdb_c7_cART_untreat.key, msigdb_c7_cART_long.key, by = "ID")
n123 <- merge(overlap_cART_c7_short_long, msigdb_c7_cART_untreat.key, by = "ID")

count_n123 <- dim(merge(overlap_cART_c7_short_long, msigdb_c7_cART_untreat.key, by = "ID"))[1]
count_n124 <- dim(merge(overlap_cART_c7_untreat_short, msigdb_c7_EC_long.key, by = "ID"))[1]
count_n134 <- dim(merge(overlap_cART_c7_untreat_long, msigdb_c7_EC_long.key, by = "ID"))[1]
count_n234 <- dim(merge(overlap_cART_c7_short_long, msigdb_c7_EC_long.key, by = "ID"))[1]

count_n12 <- dim(merge(msigdb_c7_cART_untreat.key, msigdb_c7_cART_short.key, by = "ID"))[1]
count_n13 <- dim(merge(msigdb_c7_cART_untreat.key, msigdb_c7_cART_long.key, by = "ID"))[1]
count_n14 <- dim(merge(msigdb_c7_cART_untreat.key, msigdb_c7_EC_long.key, by = "ID"))[1]
count_n23 <- dim(merge(msigdb_c7_cART_short.key, msigdb_c7_cART_long.key, by = "ID"))[1]
count_n24 <- dim(merge(msigdb_c7_cART_short.key, msigdb_c7_EC_long.key, by = "ID"))[1]
count_n34 <- dim(merge(msigdb_c7_cART_long.key, msigdb_c7_EC_long.key, by = "ID"))[1]

count_n1234 <- dim(merge(n123, msigdb_c7_EC_long.key, by = "ID"))[1]

pdf("/media/chen/DATA/evoPath/Abb/msigdb_cART_EC_c7_venn_diagram.pdf")  
draw.quad.venn(dim(msigdb_c7_cART_untreat.key)[1], dim(msigdb_c7_cART_short.key)[1], dim(msigdb_c7_cART_long.key)[1], dim(msigdb_c7_EC_long.key)[1], count_n12, count_n13, count_n14, count_n23, count_n24, count_n34, count_n123, count_n124, count_n134, count_n234, count_n1234, category = c("cART untreat", "cART short", "cART long", "EC long"))
dev.off()

# 3. r unique C7 GS
msigdb_c7_cART_short_uniq <- dplyr::anti_join(msigdb_c7_cART_short.key, msigdb_c7_cART_long.key, by = "ID")
msigdb_c7_cART_short_uniq <- dplyr::anti_join(msigdb_c7_cART_short_uniq, msigdb_c7_cART_untreat.key, by = "ID")
write.table(msigdb_c7_cART_short_uniq, file = "/media/chen/DATA/evoPath/df/C7/msigdb_c7_cART_short_uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)

msigdb_c7_cART_long_uniq <- dplyr::anti_join(msigdb_c7_cART_long.key, msigdb_c7_cART_short.key, by = "ID")
msigdb_c7_cART_long_uniq <- dplyr::anti_join(msigdb_c7_cART_long_uniq, msigdb_c7_cART_untreat.key, by = "ID")
write.table(msigdb_c7_cART_long_uniq, file = "/media/chen/DATA/evoPath/df/C7/msigdb_c7_cART_long_uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)

***include EC***
### include EC
c7_cART_untreat_uniq <- dplyr::anti_join(msigdb_c7_cART_untreat.key, msigdb_c7_cART_short.key, by = "ID")
c7_cART_untreat_uniq <- dplyr::anti_join(c7_cART_untreat_uniq, msigdb_c7_cART_long.key, by = "ID")
c7_cART_untreat_uniq <- dplyr::anti_join(c7_cART_untreat_uniq, msigdb_c7_EC_long.key, by = "ID")
write.table(c7_cART_untreat_uniq, file = "/media/chen/DATA/evoPath/df/C7/c7_cART_untreat_uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)

c7_cART_short_uniq <- dplyr::anti_join(msigdb_c7_cART_short.key, msigdb_c7_cART_untreat.key, by = "ID")
c7_cART_short_uniq <- dplyr::anti_join(c7_cART_short_uniq, msigdb_c7_cART_long.key, by = "ID")
c7_cART_short_uniq <- dplyr::anti_join(c7_cART_short_uniq, msigdb_c7_EC_long.key, by = "ID")
write.table(c7_cART_short_uniq, file = "/media/chen/DATA/evoPath/df/C7/c7_cART_short_uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)

c7_cART_long_uniq <- dplyr::anti_join(msigdb_c7_cART_long.key, msigdb_c7_cART_untreat.key, by = "ID")
c7_cART_long_uniq <- dplyr::anti_join(c7_cART_long_uniq, msigdb_c7_cART_short.key, by = "ID")
c7_cART_long_uniq <- dplyr::anti_join(c7_cART_long_uniq, msigdb_c7_EC_long.key, by = "ID")
write.table(c7_cART_long_uniq, file = "/media/chen/DATA/evoPath/df/C7/c7_cART_long_uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)

c7_EC_long_uniq <- dplyr::anti_join(msigdb_c7_EC_long.key, msigdb_c7_cART_untreat.key, by = "ID")
c7_EC_long_uniq <- dplyr::anti_join(c7_EC_long_uniq, msigdb_c7_cART_short.key, by = "ID")
c7_EC_long_uniq <- dplyr::anti_join(c7_EC_long_uniq, msigdb_c7_cART_long.key, by = "ID")
write.table(c7_EC_long_uniq, file = "/media/chen/DATA/evoPath/df/C7/c7_EC_long_uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# 4. categorize the unique paths from cARAT short & cART long
category <- c("CD4 T cells", "CD8 T cells", "B cells", "DC/macrophages/monocytes", "plasma cells", "NK cells", "cytokines", "chemokines", "others")
percentage <- c(9/42, 9/42, 6/42, 9/42, 1/42, 0/42, 5/42, 0/42, 11/42, 26/130, 11/130, 10/130, 27/130, 1/130, 4/130, 17/130, 3/130, 45/130)
condition <- c(rep("cART short", 9), rep("cART long", 9))

df_category_uniq_cART_short_long <- data.frame(category, percentage, condition)

deposer_cond <- c("cART short", "cART long")
df_category_uniq_cART_short_long$condition <- factor(df_category_uniq_cART_short_long$condition, levels = deposer_cond)

pdf("/media/chen/DATA/evoPath/Abb/df_category_C7_uniq_cART_short_long.bar.pdf", height = 5, width = 3.5)  
ggplot(df_category_uniq_cART_short_long, aes(x = condition, y = percentage, fill = category))+geom_bar(stat = "identity", position = "stack", line = "black")+scale_fill_brewer(type = "seq", palette = "RdYlBu")+theme_bw()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, vjust=1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

# 5. Clustering the unique paths from cARAT short & cART long
category <- c("CD4 T cells", "CD8 T cells", "B cells", "DC/macrophages/monocytes", "plasma cells", "NK cells", "cytokines", "chemokines", "others")
cART_uniq_C7GS_short <- c(9/42, 9/42, 6/42, 9/42, 1/42, 0/42, 5/42, 0/42, 11/42)
cART_uniq_C7GS_long <- c(26/130, 11/130, 10/130, 27/130, 1/130, 4/130, 17/130, 3/130, 45/130)
df_category_uniq_cART_short_long_pair <- data.frame(category, cART_uniq_C7GS_short, cART_uniq_C7GS_long)
row.names(df_category_uniq_cART_short_long_pair) <- df_category_uniq_cART_short_long_pair$category
df_category_uniq_cART_short_long_pair$category <- NULL
df_category_uniq_cART_short_long_pair.mx <- data.matrix(df_category_uniq_cART_short_long_pair, rownames.force = T)

# correlation plot
df_category_uniq_cART_short_long_pair.cond.co <- round(cor(df_category_uniq_cART_short_long_pair.mx), 1)
pdf("/media/chen/DATA/evoPath/Abb/df_category_uniq_cART_short_long_pair.cond.co.pdf")  
ggcorrplot(df_category_uniq_cART_short_long_pair.cond.co, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white","springgreen3"), 
           ggtheme=theme_bw)
dev.off()

pdf("/media/chen/DATA/evoPath/Abb/df_category_uniq_cART_short_long_pair.cells.co.pdf")
df_category_uniq_cART_short_long_pair.cells.co <- round(cor(t(df_category_uniq_cART_short_long_pair.mx)), 1)
ggcorrplot(scale(df_category_uniq_cART_short_long_pair.cells.co), hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white","springgreen3"), 
           ggtheme=theme_bw)
dev.off()

#PCA
df_category_uniq_cART_short_long_pair.mx.pca <- PCA(df_category_uniq_cART_short_long_pair.mx, graph = F)
eig.val <- get_eigenvalue(df_category_uniq_cART_short_long_pair.mx.pca)
pdf("/media/chen/DATA/evoPath/Abb/df_category_uniq_cART_short_long_pair.cells.pca.pdf")
fviz_pca_ind(df_category_uniq_cART_short_long_pair.mx.pca, pointsize = "cos2", pointshape = 21, fill = "#E7B800", repel = T)
dev.off()

# Slope plot
cART_uniq_C7GS_short <- c(9/42, 9/42, 6/42, 9/42, 1/42, 0/42, 5/42, 0/42, 11/42)
cART_uniq_C7GS_long <- c(26/130, 11/130, 10/130, 27/130, 1/130, 4/130, 17/130, 3/130, 45/130)
df_category_uniq_cART_short_long_pair <- data.frame(category, cART_uniq_C7GS_short, cART_uniq_C7GS_long)

left_label <- paste(df_category_uniq_cART_short_long_pair$category)
right_label <- paste(df_category_uniq_cART_short_long_pair$category)

pdf("/media/chen/DATA/evoPath/Abb/df_category_uniq_cART_short_long_pair.cells.slop.pdf")
ggplot(df_category_uniq_cART_short_long_pair, aes(x = 1, xend = 2 , y = cART_uniq_C7GS_short, yend = cART_uniq_C7GS_long, col = category), size = .75, show.legend = F)+geom_segment()+theme_bw()+geom_vline(xintercept=1, linetype="dashed", size=.1)+geom_vline(xintercept=2, linetype="dashed", size=.1)+scale_color_manual(labels = c("Up", "Down"), values = c("green"="#00ba38", "red"="#f8766d"))+labs(x = "", y = "Frequency of terms found in pathways")+xlim(.5, 2.5) + ylim(0,(1.1*(max(df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_short, df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_long))))+geom_text(label=left_label, y=df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_short, x=rep(1, NROW(df_category_uniq_cART_short_long_pair)), hjust=1.1, size=3.5)+geom_text(label=right_label, y=df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_long, x=rep(2, NROW(df_category_uniq_cART_short_long_pair)), hjust=-0.1, size=3.5)+geom_text(label="cART short", x=1, y=1.1*(max(df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_short, df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_long)), hjust=1.2, size=5)+geom_text(label="cART long", x=2, y=1.1*(max(df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_short, df_category_uniq_cART_short_long_pair$cART_uniq_C7GS_long)), hjust=-0.1, size=5)+theme(panel.background = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), panel.border = element_blank(), plot.margin = unit(c(1,2,1,2), "cm"), legend.position = "none")
dev.off()

# 6. r clustering cytokines
cART_EC_cytokines_used <- read.table("/media/chen/DATA/evoPath/df/C7/c7_cART_EC_pt_cytokines.tsv", header = T, stringsAsFactors = F)
rownames(cART_EC_cytokines_used) <- cART_EC_cytokines_used$Cytokine
cART_EC_cytokines_used$Cytokine <- NULL
cART_EC_cytokines_used.mx <- data.matrix(cART_EC_cytokines_used, rownames.force = T)

# heatmap
pdf("/media/chen/DATA/evoPath/Abb/msigdb_cART_EC_pt_cytokines_used.pdf", height = 5, width = 4)  
Heatmap(cART_EC_cytokines_used.mx, name = "Frequency", rect_gp = gpar(col = "white", lwd = 0.5), column_title = "clinical condition", row_title = "Cytokines", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), row_km = 4)
dev.off()
