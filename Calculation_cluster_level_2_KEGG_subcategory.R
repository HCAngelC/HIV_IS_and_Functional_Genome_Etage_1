```{r calculate the percentage of KEGG subcat per signal}
calculation_KEGG_subcategory <- function(df_short, df_long) {
  df_short <- dplyr::select(df_short, KEGG_subcat, count)
  df_short.agg <- aggregate(count ~ KEGG_subcat, data = df_short, FUN = sum) %>% dplyr::mutate(percent = (count/dim(df_short)[1])*100, period = "short")

  df_long <- dplyr::select(df_long, KEGG_subcat, count)
  df_long.agg <- aggregate(count ~ KEGG_subcat, data = df_long, FUN = sum)
  df_long.agg <- aggregate(count ~ KEGG_subcat, data = df_long, FUN = sum) %>% dplyr::mutate(percent = (count/dim(df_long)[1])*100, period = "long")
  
  df_tous <- dplyr::bind_rows(df_short.agg, df_long.agg)
}

calcu_KEGG_subcat_TGFB <- calculation_KEGG_subcategory(c7_kegg_short_TGFB, c7_kegg_long_TGFB) %>% dplyr::mutate(signal = "TGFB")
calcu_KEGG_subcat_IFNB <- calculation_KEGG_subcategory(c7_kegg_short_IFNB, c7_kegg_long_IFNB) %>% dplyr::mutate(signal = "IFNB")
calcu_KEGG_subcat_IFNG <- calculation_KEGG_subcategory(c7_kegg_short_IFNG, c7_kegg_long_IFNG) %>% dplyr::mutate(signal = "IFNG")
calcu_KEGG_subcat_IL10 <- calculation_KEGG_subcategory(c7_kegg_short_IL10, c7_kegg_long_IL10) %>% dplyr::mutate(signal = "IL10")
calcu_KEGG_subcat_FOXP3 <- calculation_KEGG_subcategory(c7_kegg_short_FOXP3, c7_kegg_long_FOXP3) %>% dplyr::mutate(signal = "FOXP3")

deposer_period <- c("short", "long")

#TGFB
colourCount_TGFB = length(unique(calcu_KEGG_subcat_TGFB$KEGG_subcat))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))

calcu_KEGG_subcat_TGFB$period <- factor(calcu_KEGG_subcat_TGFB$period, levels = deposer_period)

pdf("/media/chen/DATA/evoPath/Abb/calcu_KEGG_subcat_TGFB.pdf", height = 5, width = 3.5)  
ggplot(calcu_KEGG_subcat_TGFB, aes(x = period, y = percent, fill = KEGG_subcat))+geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_manual(values = getPalette(colourCount_TGFB))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

#IFNB
colourCount_IFNB = length(unique(calcu_KEGG_subcat_IFNB$KEGG_subcat))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))

calcu_KEGG_subcat_IFNB$period <- factor(calcu_KEGG_subcat_IFNB$period, levels = deposer_period)

pdf("/media/chen/DATA/evoPath/Abb/calcu_KEGG_subcat_IFNB.pdf", height = 5, width = 3.5)  
ggplot(calcu_KEGG_subcat_IFNB, aes(x = period, y = percent, fill = KEGG_subcat))+geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_manual(values = getPalette(colourCount_IFNB))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

#IFNG
colourCount_IFNG = length(unique(calcu_KEGG_subcat_IFNG$KEGG_subcat))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))

calcu_KEGG_subcat_IFNG$period <- factor(calcu_KEGG_subcat_IFNG$period, levels = deposer_period)

pdf("/media/chen/DATA/evoPath/Abb/calcu_KEGG_subcat_IFNG.pdf", height = 5, width = 4)  
ggplot(calcu_KEGG_subcat_IFNG, aes(x = period, y = percent, fill = KEGG_subcat))+geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_manual(values = getPalette(colourCount_IFNG))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

#IL10
colourCount_IL10 = length(unique(calcu_KEGG_subcat_IL10$KEGG_subcat))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))

calcu_KEGG_subcat_IL10$period <- factor(calcu_KEGG_subcat_IL10$period, levels = deposer_period)

pdf("/media/chen/DATA/evoPath/Abb/calcu_KEGG_subcat_IL10.pdf", height = 5, width = 4)  
ggplot(calcu_KEGG_subcat_IL10, aes(x = period, y = percent, fill = KEGG_subcat))+geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_manual(values = getPalette(colourCount_IL10))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()

#FOXP3
colourCount_FOXP3 = length(unique(calcu_KEGG_subcat_FOXP3$KEGG_subcat))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))

calcu_KEGG_subcat_FOXP3$period <- factor(calcu_KEGG_subcat_FOXP3$period, levels = deposer_period)

pdf("/media/chen/DATA/evoPath/Abb/calcu_KEGG_subcat_FOXP3.pdf", height = 5, width = 6.5)  
ggplot(calcu_KEGG_subcat_FOXP3, aes(x = period, y = percent, fill = KEGG_subcat))+geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_manual(values = getPalette(colourCount_FOXP3))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()
