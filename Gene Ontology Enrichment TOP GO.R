# Gene Ontology Enrichment Analysis using topGO

# 0. Install and Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("org.Hs.eg.db", "topGO", "ggplot2", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(org.Hs.eg.db)
library(topGO)
library(ggplot2)
library(ggrepel)

# 1. Load your DEG table (replace with your actual file)
deg <- read.csv("C:/Users/kanis/Downloads/Differential Expression Analysis - EdgeR/DE_results_edgeR.csv")
head(deg)

# 2. Map your IDs to gene symbols
entrez_ids <- as.character(deg$X)
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = entrez_ids,
                                      column = "SYMBOL",
                                      keytype = "ENTREZID",
                                      multiVals = "first")

deg$GeneSymbol <- gene_symbols
deg <- deg[!is.na(deg$GeneSymbol), ]
head(deg)

# 3. Save annotated table (optional)
write.csv(deg, "DEGs_with_GeneSymbols.csv", row.names = FALSE)

# 4. Prepare gene list for topGO
geneList <- deg$FDR     # or deg$PValue if FDR is not available
names(geneList) <- deg$GeneSymbol

# 5. Define function for significant genes
topDiffGenes <- function(pval) { return(pval < 0.05) }

# 6. topGO enrichment for BP, CC, MF
GOdata_BP <- new("topGOdata",
                 description = "GO enrichment - BP",
                 ontology = "BP",
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "symbol",
                 nodeSize = 10)
GOdata_CC <- new("topGOdata",
                 description = "GO enrichment - CC",
                 ontology = "CC",
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "symbol",
                 nodeSize = 10)
GOdata_MF <- new("topGOdata",
                 description = "GO enrichment - MF",
                 ontology = "MF",
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org,
                 mapping = "org.Hs.eg.db",
                 ID = "symbol",
                 nodeSize = 10)

# 7. Run Fisher's exact test for each ontology
resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")
resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")

# 8. Extract top results
run_topGO <- function(GOdata, resultFisher, num = 10, ontology = "BP") {
  tab <- GenTable(GOdata,
                  classicFisher = resultFisher,
                  orderBy = "classicFisher",
                  topNodes = num)
  tab$Ontology <- ontology
  tab$classicFisher <- as.numeric(tab$classicFisher)
  tab
}
res_BP <- run_topGO(GOdata_BP, resultFisher_BP, 10, "BP")
res_CC <- run_topGO(GOdata_CC, resultFisher_CC, 10, "CC")
res_MF <- run_topGO(GOdata_MF, resultFisher_MF, 10, "MF")

go_all <- rbind(res_BP, res_CC, res_MF)
go_all <- go_all[!is.na(go_all$classicFisher) & go_all$classicFisher > 0, ]

# 9. Plotting: Bar plot
bar_plot <- ggplot(go_all, aes(x = reorder(Term, -log10(classicFisher)),
                               y = -log10(classicFisher),
                               fill = Ontology)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  coord_flip() +
  labs(title = "Top Enriched GO Terms (BP, CC, MF)",
       x = "GO Term",
       y = expression(-log[10]*"(p-value)")) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("BP" = "#1b9e77", "CC" = "#7570b3", "MF" = "#d95f02")) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
print(bar_plot)
ggsave("GO_enrichment_barplot.png", bar_plot, width = 10, height = 6, dpi = 300)

# 10. Plotting: Bubble plot
bubble_plot <- ggplot(go_all, aes(x = Ontology,
                                  y = reorder(Term, -log10(classicFisher)),
                                  size = Significant,
                                  color = -log10(classicFisher))) +
  geom_point(alpha = 0.7) +
  ggrepel::geom_text_repel(aes(label = Term), size = 3, max.overlaps = 10) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 10)) +
  labs(title = "GO Bubble Plot (BP, CC, MF)",
       x = "Ontology",
       y = "GO Term",
       size = "Significant Genes",
       color = "-log10(p-value)") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
print(bubble_plot)
ggsave("GO_bubbleplot_white_bg.png", bubble_plot, width = 12, height = 8, dpi = 300, bg = "white")
