### 1. import IS von Einkauf et al. 2022###
cART_untreat <- read.table("/media/chen/DATA/evoPath/df/IS_cART_Einkauf_et_al_2022/IS_untreat_cART.csv", header = F, stringsAsFactors = F) %>% unique()
cART_short <- read.table("/media/chen/DATA/evoPath/df/IS_cART_Einkauf_et_al_2022/IS_short_cART.csv", header = F, stringsAsFactors = F) %>% unique()
cART_long <- read.table("/media/chen/DATA/evoPath/df/IS_cART_Einkauf_et_al_2022/IS_long_cART.csv", header = F, stringsAsFactors = F) %>% unique()

### 2. obtenir ENTREZID ###
cART_untreat$V1 <- as.character(cART_untreat$V1)
cART_untreat.enID <- bitr(cART_untreat$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

cART_short$V1 <- as.character(cART_short$V1)
cART_short.enID <- bitr(cART_short$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

cART_long$V1 <- as.character(cART_long$V1)
cART_long.enID <- bitr(cART_long$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

### 3. faire des lists des ENTRZID & comparaison GO & KEGG ###
List_tous_cART <- list(untreat_cART = cART_untreat.enID$ENTREZID,  short_cART = cART_short.enID$ENTREZID, long_cART = cART_long.enID$ENTREZID)

#Comp_tous_cART_GO <- compareCluster(geneClusters = List_tous_cART, fun = enrichGO)
Comp_tous_cART_KEGG <- compareCluster(geneClusters = List_tous_cART, fun = enrichKEGG)
Comp_tous_cART_KEGG <- setReadable(Comp_tous_cART_KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

### 4. plot
pdf("/media/chen/DATA/evoPath/Abb/Enriched_Pw_tous_cART_Einkauf_et_al_2022.pdf", height = 5, width = 5.5)
dotplot(Comp_tous_cART_KEGG)
dev.off()

### 5. enrichKEGG
enkegg_cART_untreat <- data.frame(enrichKEGG(gene = cART_untreat.enID$ENTREZID, organism = "hsa", pvalueCutoff = 0.5))
enkegg_cART_short <- data.frame(enrichKEGG(gene = cART_short.enID$ENTREZID, organism = "hsa", pvalueCutoff = 0.5))
enkegg_cART_long <- data.frame(enrichKEGG(gene = cART_long.enID$ENTREZID, organism = "hsa", pvalueCutoff = 0.5))

### 6. clustering enriched KEGG
enkegg_cART_untreat_key <- enkegg_cART_untreat %>% dplyr::select(ID, Description, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "untreat")
enkegg_cART_untreat_key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", enkegg_cART_untreat_key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", enkegg_cART_untreat_key$GeneRatio, perl=T))
enkegg_cART_untreat_key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", enkegg_cART_untreat_key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", enkegg_cART_untreat_key$BgRatio, perl=T))
## convert fraction GeneRatio & BgRatio to numeric by using the pachage "gsubfn" and calculate the ratio <- GeneRatio/BgRatio

enkegg_cART_long_key <- enkegg_cART_long %>% dplyr::select(ID, Description, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "long")
enkegg_cART_long_key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", enkegg_cART_long_key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", enkegg_cART_long_key$GeneRatio, perl=T))
enkegg_cART_long_key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", enkegg_cART_long_key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", enkegg_cART_long_key$BgRatio, perl=T))

enkegg_cART_short_key <- enkegg_cART_short %>% dplyr::select(ID, Description, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "short")
enkegg_cART_short_key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", enkegg_cART_short_key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", enkegg_cART_short_key$GeneRatio, perl=T))
enkegg_cART_short_key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", enkegg_cART_short_key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", enkegg_cART_short_key$BgRatio, perl=T))

enkegg_cART_tous <- dplyr::bind_rows(enkegg_cART_untreat_key, enkegg_cART_long_key, enkegg_cART_short_key) %>% dplyr::mutate(Ratio = GeneRatio/BgRatio)
enkegg_cART_tous.untreat <- dplyr::filter(enkegg_cART_tous, cond == "untreat")
enkegg_cART_tous.long <- dplyr::filter(enkegg_cART_tous, cond == "long")
enkegg_cART_tous.short <- dplyr::filter(enkegg_cART_tous, cond == "short")

enkegg_cART_tous_ID_unique <- enkegg_cART_tous %>% dplyr::select(ID) %>% unique()

enkegg_cART_tous.untreat_avoir <- dplyr::inner_join(enkegg_cART_tous.untreat, enkegg_cART_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
enkegg_cART_tous.untreat_sans <- dplyr::anti_join(enkegg_cART_tous_ID_unique, enkegg_cART_tous.untreat, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
enkegg_cART_tous.untreat.ID.Ratio <- dplyr::bind_rows(enkegg_cART_tous.untreat_avoir, enkegg_cART_tous.untreat_sans) %>% dplyr::rename(untreat = Ratio)

enkegg_cART_tous.long_avoir <- dplyr::inner_join(enkegg_cART_tous.long, enkegg_cART_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
enkegg_cART_tous.long_sans <- dplyr::anti_join(enkegg_cART_tous_ID_unique, enkegg_cART_tous.long, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
enkegg_cART_tous.long.ID.Ratio <- dplyr::bind_rows(enkegg_cART_tous.long_avoir, enkegg_cART_tous.long_sans) %>% dplyr::rename(long = Ratio)

enkegg_cART_tous.short_avoir <- dplyr::inner_join(enkegg_cART_tous.short, enkegg_cART_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
enkegg_cART_tous.short_sans <- dplyr::anti_join(enkegg_cART_tous_ID_unique, enkegg_cART_tous.short, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
enkegg_cART_tous.short.ID.Ratio <- dplyr::bind_rows(enkegg_cART_tous.short_avoir, enkegg_cART_tous.short_sans) %>% dplyr::rename(short = Ratio)

enkegg_cART_tous.ID.Ratio <- dplyr::bind_cols(enkegg_cART_tous.untreat.ID.Ratio, enkegg_cART_tous.long.ID.Ratio, enkegg_cART_tous.short.ID.Ratio) %>% dplyr::select(ID...1, untreat, short, long)
row.names(enkegg_cART_tous.ID.Ratio) <- enkegg_cART_tous.ID.Ratio$ID...1
enkegg_cART_tous.ID.Ratio$ID...1 <- NULL
enkegg_cART_tous.ID.Ratio.mx <- data.matrix(enkegg_cART_tous.ID.Ratio, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/enkegg_cART_evolution_pathway.pdf")  
Heatmap(enkegg_cART_tous.ID.Ratio.mx, name = "GeneRatio/BgRatio", rect_gp = gpar(col = "white", lwd = 0.5), column_title = "clinical condition", row_title = "KEGG pathways", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), row_km = 4)
dev.off()
