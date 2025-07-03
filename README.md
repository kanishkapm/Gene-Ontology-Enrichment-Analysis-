# Gene-Ontology-Enrichment-Analysis-

R workflow for DEG annotation and GO enrichment with visualization. Scripts for differential expression and GO analysis in R. R analysis pipeline for DEGs and Gene Ontology enrichment. Automated DEG and GO term analysis with R, topGO, and clusterProfiler. Tools for RNA-seq DEG annotation, GO enrichment, and visualization in R.

# Human Differential Expression and GO Enrichment Workflow

This repository contains an R workflow for performing differential expression analysis and GO enrichment (BP, CC, MF) with visualization.

## Requirements

- R (>= 4.0)
- Bioconductor packages: `org.Hs.eg.db`, `topGO`, `clusterProfiler`
- CRAN packages: `ggplot2`, `ggrepel`

## Usage

1. Edit the file paths in `DEA_GO_enrichment.R` to point to your own DEG results file.
2. Run the script in R or RStudio.
3. The script will produce annotated DEGs, GO enrichment results (using both topGO and clusterProfiler), and visualization plots.

## Outputs

- `DEGs_with_GeneSymbols.csv` — DEGs annotated with gene symbols
- `GO_enrichment_barplot.png` — Bar plot of top GO terms (topGO)
- `GO_bubbleplot_white_bg.png` — Bubble plot of GO terms (topGO)
- `clusterProfiler_GO_enrichment.csv` — GO enrichment results table from clusterProfiler
- `clusterProfiler_GO_barplot.png` — Bar plot of GO terms from clusterProfiler
- `clusterProfiler_GO_dotplot.png` — Dot plot of GO terms from clusterProfiler

## License

[MIT](LICENSE)
