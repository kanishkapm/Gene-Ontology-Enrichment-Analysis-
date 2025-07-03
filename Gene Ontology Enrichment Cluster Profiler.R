# Gene Ontology Enrichment Analysis using clusterProfiler

# 0. Install and Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 1. Load your DEG table (replace with your actual file)
deg <- read.csv("C:/Users/kanis/Downloads/Differential Expression Analysis - EdgeR/DE_results_edgeR.csv")
entrez_ids <- as.character(deg$X)

# 2. Run GO enrichment for BP, MF, CC
ego_bp <- enrichGO(gene         = entrez_ids,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "BP",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable     = TRUE)

ego_mf <- enrichGO(gene         = entrez_ids,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "MF",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable     = TRUE)

ego_cc <- enrichGO(gene         = entrez_ids,
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   ont          = "CC",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable     = TRUE)

# 3. Visualize results (dotplot)
dotplot(ego_bp, showCategory=10, title="GO Enrichment: Biological Process")
dotplot(ego_mf, showCategory=10, title="GO Enrichment: Molecular Function")
dotplot(ego_cc, showCategory=10, title="GO Enrichment: Cellular Component")

# Save result tables
write.csv(as.data.frame(ego_bp), "GO_BP_enrichment_clusterProfiler.csv", row.names = FALSE)
write.csv(as.data.frame(ego_mf), "GO_MF_enrichment_clusterProfiler.csv", row.names = FALSE)
write.csv(as.data.frame(ego_cc), "GO_CC_enrichment_clusterProfiler.csv", row.names = FALSE)
