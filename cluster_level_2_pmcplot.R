kegg_cART_short_TGFB_rel_genes_a <- enrichKEGG(gene = cART_short_TGFB$V1, organism = "hsa", pvalueCutoff = 0.2)
kegg_cART_long_TGFB_rel_genes_a <- enrichKEGG(gene = cART_long_TGFB$V1, organism = "hsa", pvalueCutoff = 0.2)

kegg_cART_short_IFNB_rel_genes_a <- enrichKEGG(gene = cART_short_IFNB$V1, organism = "hsa", pvalueCutoff = 0.2)
kegg_cART_long_IFNB_rel_genes_a <- enrichKEGG(gene = cART_long_IFNB$V1, organism = "hsa", pvalueCutoff = 0.2)

kegg_cART_short_IFNG_rel_genes_a <- enrichKEGG(gene = cART_short_IFNG$V1, organism = "hsa", pvalueCutoff = 0.2)
kegg_cART_long_IFNG_rel_genes_a <- enrichKEGG(gene = cART_long_IFNG$V1, organism = "hsa", pvalueCutoff = 0.2)

kegg_cART_short_IL10_rel_genes_a <- enrichKEGG(gene = cART_short_IL10$V1, organism = "hsa", pvalueCutoff = 0.2)
kegg_cART_long_IL10_rel_genes_a <- enrichKEGG(gene = cART_long_IL10$V1, organism = "hsa", pvalueCutoff = 0.2)

kegg_cART_short_FOXP3_rel_genes_a <- enrichKEGG(gene = cART_short_FOXP3$V1, organism = "hsa", pvalueCutoff = 0.2)
kegg_cART_long_FOXP3_rel_genes_a <- enrichKEGG(gene = cART_long_FOXP3$V1, organism = "hsa", pvalueCutoff = 0.2)

# Terms used in PMC
terms_kegg_cART_short_TGFB <- kegg_cART_short_TGFB_rel_genes_a$Description
terms_kegg_cART_short_TGFB_plot <- enrichplot::pmcplot(terms_kegg_cART_short_TGFB, 2010:2021, proportion=FALSE)
terms_kegg_cART_long_TGFB <- kegg_cART_long_TGFB_rel_genes_a$Description
terms_kegg_cART_long_TGFB_plot <- enrichplot::pmcplot(terms_kegg_cART_long_TGFB, 2010:2021, proportion=FALSE)

terms_kegg_cART_short_IFNB <- kegg_cART_short_IFNB_rel_genes_a$Description
terms_kegg_cART_short_IFNB_plot <- enrichplot::pmcplot(terms_kegg_cART_short_IFNB, 2010:2021, proportion=FALSE)
terms_kegg_cART_long_IFNB <- kegg_cART_long_IFNB_rel_genes_a$Description
terms_kegg_cART_long_IFNB_plot <- enrichplot::pmcplot(terms_kegg_cART_long_IFNB, 2010:2021, proportion=FALSE)

terms_kegg_cART_short_IFNG <- kegg_cART_short_IFNG_rel_genes_a$Description
terms_kegg_cART_short_IFNG_plot <- enrichplot::pmcplot(terms_kegg_cART_short_IFNG, 2010:2021, proportion=FALSE)
terms_kegg_cART_long_IFNG <- kegg_cART_long_IFNG_rel_genes_a$Description
terms_kegg_cART_long_IFNG_plot <- enrichplot::pmcplot(terms_kegg_cART_long_IFNG, 2010:2021, proportion=FALSE)

terms_kegg_cART_short_IL10 <- kegg_cART_short_IL10_rel_genes_a$Description
terms_kegg_cART_short_IL10_plot <- enrichplot::pmcplot(terms_kegg_cART_short_IL10, 2010:2021, proportion=FALSE)
terms_kegg_cART_long_IL10 <- kegg_cART_long_IL10_rel_genes_a$Description
terms_kegg_cART_long_IL10_plot <- enrichplot::pmcplot(terms_kegg_cART_long_IL10, 2010:2021, proportion=FALSE)

terms_kegg_cART_short_FOXP3 <- kegg_cART_short_FOXP3_rel_genes_a$Description
terms_kegg_cART_short_FOXP3_plot <- enrichplot::pmcplot(terms_kegg_cART_short_FOXP3, 2010:2021, proportion=FALSE)
terms_kegg_cART_long_FOXP3 <- kegg_cART_long_FOXP3_rel_genes_a$Description
terms_kegg_cART_long_FOXP3_plot <- enrichplot::pmcplot(terms_kegg_cART_long_FOXP3, 2010:2021, proportion=FALSE)
