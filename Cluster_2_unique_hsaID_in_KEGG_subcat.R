### *** sort out unique hasID for short- & long cART*** ###
c7_kegg_short_subcat_term_tous <- dplyr::bind_rows(c7_kegg_short_TGFB, c7_kegg_short_IFNB, c7_kegg_short_IFNG, c7_kegg_short_IL10, c7_kegg_short_CXCR5) %>% dplyr::select(ID, KEGG_subcat) %>% unique()
c7_kegg_long_subcat_term_tous <- dplyr::bind_rows(c7_kegg_long_TGFB, c7_kegg_long_IFNB, c7_kegg_long_IFNG, c7_kegg_long_IL10, c7_kegg_long_CXCR5) %>% dplyr::select(ID, KEGG_subcat) %>% unique()

c7_kegg_short_ID_sans_FOXP3.agg.uniq <- dplyr::anti_join(c7_kegg_short_ID_sans_FOXP3.agg, c7_kegg_long_ID_sans_FOXP3.agg, by = "ID")
c7_kegg_long_ID_sans_FOXP3.agg.uniq <- dplyr::anti_join(c7_kegg_long_ID_sans_FOXP3.agg, c7_kegg_short_ID_sans_FOXP3.agg, by = "ID")

c7_kegg_short_subcat_uniq_tous <- dplyr::inner_join(c7_kegg_short_ID_sans_FOXP3.agg.uniq, c7_kegg_short_subcat_term_tous, by = "ID")
c7_kegg_long_subcat_uniq_tous <- dplyr::inner_join(c7_kegg_long_ID_sans_FOXP3.agg.uniq, c7_kegg_long_subcat_term_tous, by = "ID")

#verify <- dplyr::inner_join(c7_kegg_short_subcat_uniq_tous, c7_kegg_long_subcat_uniq_tous, by = "ID")

pdf("/media/chen/DATA/evoPath/Abb/c7_kegg_short_subcat_uniq_tous.pdf", height = 6, width = 7)
ggplot(c7_kegg_short_subcat_uniq_tous, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(44))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/c7_kegg_short_subcat_uniq_tous.svg", height = 6, width = 7)
ggplot(c7_kegg_short_subcat_uniq_tous, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Blues"))(44))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

pdf("/media/chen/DATA/evoPath/Abb/c7_kegg_long_subcat_uniq_tous.pdf", height = 6, width = 6.5)
ggplot(c7_kegg_long_subcat_uniq_tous, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Reds"))(25))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

svg("/media/chen/DATA/evoPath/Abb/c7_kegg_long_subcat_uniq_tous.svg", height = 6, width = 6.5)
ggplot(c7_kegg_long_subcat_uniq_tous, aes(x = KEGG_subcat, y = count, fill = ID))+geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Reds"))(25))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()
