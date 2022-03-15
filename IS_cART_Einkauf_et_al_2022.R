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
