### 1. import IS ###
cART_untreat <- read.table("/media/chen/LaCie/IGH_backup/evoPath/input_genelist/IS_untreat_cART.csv", header = F, stringsAsFactors = F) %>% unique()
cART_short <- read.table("/media/chen/LaCie/IGH_backup/evoPath/input_genelist/IS_short_cART.csv", header = F, stringsAsFactors = F) %>% unique()
cART_long <- read.table("/media/chen/LaCie/IGH_backup/evoPath/input_genelist/IS_long_cART.csv", header = F, stringsAsFactors = F) %>% unique()
long_term_EC <- read.table("/media/chen/LaCie/IGH_backup/evoPath/input_genelist/IS_Jiang_et_al_2021_Gene_list.txt", header = F, stringsAsFactors = F) %>% dplyr::select(V1) %>% unique()

### 2. obtenir ENTREZID ###
cART_untreat$V1 <- as.character(cART_untreat$V1)
cART_untreat.enID <- bitr(cART_untreat$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

cART_short$V1 <- as.character(cART_short$V1)
cART_short.enID <- bitr(cART_short$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

cART_long$V1 <- as.character(cART_long$V1)
cART_long.enID <- bitr(cART_long$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

long_term_EC$V1 <- as.character(long_term_EC$V1)
long_term_EC.enID <- bitr(long_term_EC$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

### 3. MSiGDb
msigdb_c7_cART_untreat <- data.frame(enricher(cART_untreat.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_cART_short <- data.frame(enricher(cART_short.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_cART_long <- data.frame(enricher(cART_long.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_c7_EC_long <- data.frame(enricher(long_term_EC.enID$ENTREZID, TERM2GENE = m_tC7))

write.table(msigdb_c7_cART_untreat, file = "/media/chen/LaCie/IGH_backup/evoPath/df/wiederholen/msigdb_c7_ART_untreat.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(msigdb_c7_cART_short, file = "/media/chen/LaCie/IGH_backup/evoPath/df/wiederholen/msigdb_c7_ART_short.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(msigdb_c7_cART_long, file = "/media/chen/LaCie/IGH_backup/evoPath/df/wiederholen/msigdb_c7_ART_long.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(msigdb_c7_EC_long, file = "/media/chen/LaCie/IGH_backup/evoPath/df/wiederholen/msigdb_c7_EC_long.txt", row.names = F, col.names = T, sep = "\t", quote = F)

msigdb_c7_cART_untreat.key <- msigdb_c7_cART_untreat %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "untreat")
msigdb_c7_cART_untreat.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_untreat.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_untreat.key$GeneRatio, perl=T))
msigdb_c7_cART_untreat.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_untreat.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_untreat.key$BgRatio, perl=T))
msigdb_c7_cART_untreat.key <- msigdb_c7_cART_untreat.key %>% dplyr::mutate(RF = GeneRatio/BgRatio)

msigdb_c7_cART_short.key <- msigdb_c7_cART_short %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "short")
msigdb_c7_cART_short.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$GeneRatio, perl=T))
msigdb_c7_cART_short.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_short.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_short.key$BgRatio, perl=T))
msigdb_c7_cART_short.key <- msigdb_c7_cART_short.key %>% dplyr::mutate(RF = GeneRatio/BgRatio)

msigdb_c7_cART_long.key <- msigdb_c7_cART_long %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "long")
msigdb_c7_cART_long.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$GeneRatio, perl=T))
msigdb_c7_cART_long.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_cART_long.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_cART_long.key$BgRatio, perl=T))
msigdb_c7_cART_long.key <- msigdb_c7_cART_long.key %>% dplyr::mutate(RF = GeneRatio/BgRatio)

msigdb_c7_EC_long.key <- msigdb_c7_EC_long %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "EC long")
msigdb_c7_EC_long.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_EC_long.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_EC_long.key$GeneRatio, perl=T))
msigdb_c7_EC_long.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_c7_EC_long.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_c7_EC_long.key$BgRatio, perl=T))
msigdb_c7_EC_long.key <- msigdb_c7_EC_long.key %>% dplyr::mutate(RF = GeneRatio/BgRatio)

df_msigdb_tous.key <- dplyr::bind_rows(msigdb_c7_cART_untreat.key, msigdb_c7_cART_short.key, msigdb_c7_cART_long.key, msigdb_c7_EC_long.key) %>% dplyr::select(cond, RF)

### 4. Random ENTREZID machen
hg38_genelist <- read.table("/media/chen/DATA/hg38_gene_list.gene_name.csv", header = F, stringsAsFactors = F) %>% unique()
hg38_genelist$V1 <- as.character(hg38_genelist$V1)
hg38_genelist.enID <- bitr(hg38_genelist$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_150_1 <- sample_n(hg38_genelist, 150, replace = F)
random_contl_150_1$V1 <- as.character(random_contl_150_1$V1)
random_contl_150_1.enID <- bitr(random_contl_150_1$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_150_2 <- sample_n(hg38_genelist, 150, replace = F)
random_contl_150_2$V1 <- as.character(random_contl_150_2$V1)
random_contl_150_2.enID <- bitr(random_contl_150_2$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_250_1 <- sample_n(hg38_genelist, 250, replace = F)
random_contl_250_1$V1 <- as.character(random_contl_250_1$V1)
random_contl_250_1.enID <- bitr(random_contl_250_1$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_250_2 <- sample_n(hg38_genelist, 250, replace = F)
random_contl_250_2$V1 <- as.character(random_contl_250_2$V1)
random_contl_250_2.enID <- bitr(random_contl_250_2$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_350 <- sample_n(hg38_genelist, 350, replace = F)
random_contl_350$V1 <- as.character(random_contl_350$V1)
random_contl_350.enID <- bitr(random_contl_350$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_450 <- sample_n(hg38_genelist, 450, replace = F)
random_contl_450$V1 <- as.character(random_contl_450$V1)
random_contl_450.enID <- bitr(random_contl_450$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_550 <- sample_n(hg38_genelist, 550, replace = F)
random_contl_550$V1 <- as.character(random_contl_550$V1)
random_contl_550.enID <- bitr(random_contl_550$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_650 <- sample_n(hg38_genelist, 650, replace = F)
random_contl_650$V1 <- as.character(random_contl_650$V1)
random_contl_650.enID <- bitr(random_contl_650$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

random_contl_1000 <- sample_n(hg38_genelist, 1000, replace = F)
random_contl_1000$V1 <- as.character(random_contl_1000$V1)
random_contl_1000.enID <- bitr(random_contl_1000$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

### 5. Run NSigDb on random control
msigdb_random_contl_150_1.enID <- data.frame(enricher(random_contl_150_1.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_random_contl_150_2.enID <- data.frame(enricher(random_contl_150_2.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_random_contl_250_1.enID <- data.frame(enricher(random_contl_250_1.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_random_contl_250_2.enID <- data.frame(enricher(random_contl_250_2.enID$ENTREZID, TERM2GENE = m_tC7))

msigdb_random_contl_350.enID <- data.frame(enricher(random_contl_350.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_random_contl_450.enID <- data.frame(enricher(random_contl_450.enID$ENTREZID, TERM2GENE = m_tC7))

msigdb_random_contl_550.enID <- data.frame(enricher(random_contl_550.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_random_contl_550.key <- msigdb_random_contl_550.enID %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "conl_550")
msigdb_random_contl_550.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_random_contl_550.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_random_contl_550.key$GeneRatio, perl=T))
msigdb_random_contl_550.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_random_contl_550.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_random_contl_550.key$BgRatio, perl=T))
msigdb_random_contl_550.key <- msigdb_random_contl_550.key %>% dplyr::mutate(RF = GeneRatio/BgRatio) %>% dplyr::select(cond, RF)

msigdb_random_contl_650.enID <- data.frame(enricher(random_contl_650.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_random_contl_1000.enID <- data.frame(enricher(random_contl_1000.enID$ENTREZID, TERM2GENE = m_tC7))

msigdb_hg38_genelist.enID <- data.frame(enricher(hg38_genelist.enID$ENTREZID, TERM2GENE = m_tC7))
msigdb_hg38_genelist.key <- msigdb_hg38_genelist.enID %>% dplyr::select(ID, GeneRatio, BgRatio) %>% dplyr::mutate(cond = "hg38")
msigdb_hg38_genelist.key$GeneRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_hg38_genelist.key$GeneRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_hg38_genelist.key$GeneRatio, perl=T))
msigdb_hg38_genelist.key$BgRatio <- as.numeric( gsub("(\\d+)/(\\d+)", "\\1", msigdb_hg38_genelist.key$BgRatio, perl=T) ) / as.numeric(gsub("(\\d+)/(\\d+)", "\\2", msigdb_hg38_genelist.key$BgRatio, perl=T))
msigdb_hg38_genelist.key <- msigdb_hg38_genelist.key %>% dplyr::mutate(RF = GeneRatio/BgRatio) %>% dplyr::select(cond, RF)

df_msigdb_tous_contl_hg38.key <- dplyr::bind_rows(df_msigdb_tous.key, msigdb_random_contl_550.key, msigdb_hg38_genelist.key)

### 6. plot
deposer_cond_<- c("untreat", "short", "long", "EC long", "conl_550", "hg38")
Vergleichung_RF_ <- list(c("untreat", "hg38"), c("short", "hg38"), c("long", "hg38"), c("EC long", "hg38"))

df_msigdb_tous_contl_hg38.key$cond <- factor(df_msigdb_tous_contl_hg38.key$cond, levels = deposer_cond_)

svg("/media/chen/LaCie/IGH_backup/evoPath/Add/df_msigdb_tous_contl_hg38.key.svg", height = 5, width = 3)
ggplot(df_msigdb_tous_contl_hg38.key, aes(x = cond, y = RF, fill = cond))+geom_boxplot()+scale_fill_brewer(palette = "Blues")+stat_compare_means(comparisons = Vergleichung_RF_, label = "p.signif")+theme_bw()+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle=45, hjust = 1, vjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
dev.off()
