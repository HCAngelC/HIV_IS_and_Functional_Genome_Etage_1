m_df <- msigdbr(species = "Homo sapiens")

m_tC7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% dplyr::select(gs_name, entrez_gene)

msigdb_c7_cART_untreat <- data.frame(enricher(cART_untreat.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_cART_short <- data.frame(enricher(cART_short.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_cART_long <- data.frame(enricher(cART_long.enID$ENTREZID, TERM2GENE = m_tC7))

msigdb_c7_cART_untreat.key <- msigdb_c7_cART_untreat %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "untreat")
msigdb_c7_cART_untreat.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_untreat.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_untreat.key$GeneRatio, perl=T))
msigdb_c7_cART_untreat.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_untreat.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_untreat.key$BgRatio, perl=T))

msigdb_c7_cART_short.key <- msigdb_c7_cART_short %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "short")
msigdb_c7_cART_short.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$GeneRatio, perl=T))
msigdb_c7_cART_short.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$BgRatio, perl=T))

msigdb_c7_cART_long.key <- msigdb_c7_cART_long %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "long")
msigdb_c7_cART_long.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$GeneRatio, perl=T))
msigdb_c7_cART_long.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$BgRatio, perl=T))

msigdb_cART_tous <- dplyr::bind_rows(msigdb_c7_cART_untreat.key, msigdb_c7_cART_long.key, msigdb_c7_cART_short.key) %>% dplyr::mutate(Ratio = GeneRatio/BgRatio)
msigdb_cART_tous.untreat <- dplyr::filter(msigdb_cART_tous, cond == "untreat")
msigdb_cART_tous.long <- dplyr::filter(msigdb_cART_tous, cond == "long")
msigdb_cART_tous.short <- dplyr::filter(msigdb_cART_tous, cond == "short")

msigdb_cART_tous_ID_unique <- msigdb_cART_tous %>% dplyr::select(ID) %>% unique()

msigdb_cART_tous.untreat_avoir <- dplyr::inner_join(msigdb_cART_tous.untreat, msigdb_cART_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_tous.untreat_sans <- dplyr::anti_join(msigdb_cART_tous_ID_unique, msigdb_cART_tous.untreat, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_tous.untreat.ID.Ratio <- dplyr::bind_rows(msigdb_cART_tous.untreat_avoir, msigdb_cART_tous.untreat_sans) %>% dplyr::rename(untreat = Ratio)

msigdb_cART_tous.long_avoir <- dplyr::inner_join(msigdb_cART_tous.long, msigdb_cART_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_tous.long_sans <- dplyr::anti_join(msigdb_cART_tous_ID_unique, msigdb_cART_tous.long, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_tous.long.ID.Ratio <- dplyr::bind_rows(msigdb_cART_tous.long_avoir, msigdb_cART_tous.long_sans) %>% dplyr::rename(long = Ratio)

msigdb_cART_tous.short_avoir <- dplyr::inner_join(msigdb_cART_tous.short, msigdb_cART_tous_ID_unique, by = "ID") %>% dplyr::select(ID, Ratio)
msigdb_cART_tous.short_sans <- dplyr::anti_join(msigdb_cART_tous_ID_unique, msigdb_cART_tous.short, by = "ID") %>% dplyr::select(ID) %>% dplyr::mutate(Ratio = 0)
msigdb_cART_tous.short.ID.Ratio <- dplyr::bind_rows(msigdb_cART_tous.short_avoir, msigdb_cART_tous.short_sans) %>% dplyr::rename(short = Ratio)

msigdb_cART_tous.ID.Ratio <- dplyr::bind_cols(msigdb_cART_tous.untreat.ID.Ratio, msigdb_cART_tous.long.ID.Ratio, msigdb_cART_tous.short.ID.Ratio) %>% dplyr::select(ID...1, untreat, short, long)
row.names(msigdb_cART_tous.ID.Ratio) <- msigdb_cART_tous.ID.Ratio$ID...1
msigdb_cART_tous.ID.Ratio$ID...1 <- NULL
msigdb_cART_tous.ID.Ratio.mx <- data.matrix(msigdb_cART_tous.ID.Ratio, rownames.force = T)

pdf("/media/chen/DATA/evoPath/Abb/msigdb_cART_evolution_pathway.pdf")  
Heatmap(msigdb_cART_tous.ID.Ratio.mx, name = "GeneRatio/BgRatio", rect_gp = gpar(col = "white", lwd = 0.5), column_title = "clinical condition", row_title = "C7: immunologic signatures", column_title_side = "bottom", row_dend_width = unit(4, "cm"), row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), row_km = 5)
dev.off()
