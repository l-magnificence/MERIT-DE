# MERIT-DE
Moderated Empirical-null Robust-Input Test for Differential Expression

## Requirements
- R ≥ 4.1
- Bioconductor dependencies `edgeR`, `limma`
- ## Install dependencies
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma"))
```
## Install MERIT-DE
```r
devtools::install_github("https://github.com/l-magnificence/MERIT-DE")
```
## Quick Start
### 1. Input format

MERIT-DE needs two things:

1. **Count matrix**: **genes (rows) × samples (columns)** of **raw integer counts**
   (not normalized, not logged). 
2. **Group**: a vector of length = number of samples with **exactly two levels**
   (e.g. `"ctrl"` / `"treat"`), aligned to the matrix columns.

> Do not pre-apply TMM/CPM/log — MERIT-DE does this internally. Pass raw counts.

---
### 2. Standard bulk RNA-seq workflow
```r
library(MERITde)
counts <- as.matrix(read.csv("counts.csv", row.names = 1))
group  <- c("ctrl","ctrl","ctrl","treat","treat","treat")

res <- merit_de(counts, group)
write.csv(res, "merit_results.csv", row.names = FALSE)
sig <- subset(res, padj < 0.05)
```

---
### 3. Adjusting for confounders (batch, age, sex, ...)

Unlike the Wilcoxon test, MERIT-DE can adjust for covariates:
```r
covar <- data.frame(batch = c("b1","b1","b2","b1","b2","b2"),
                    age   = c(45,52,39,60,48,55))
res <- merit_de(counts, group, covariates = covar)
```
### 4. Single-cell pseudobulk
Do **not** run DE per cell — that inflates false positives. Aggregate (sum) cells by
individual × cell-type × condition into a pseudobulk matrix (sample unit = individual),
then run MERIT-DE **per cell type**.
```r
library(MERITde)
library(Matrix)

# counts : genes x cells raw counts (matrix or dgCMatrix)
# meta   : data.frame aligned to the COLUMNS of `counts`, with columns
#          individual, cell_type, condition (exactly two levels)

results <- list()
for (ct in unique(meta$cell_type)) {
  keep <- meta$cell_type == ct
  cnt  <- counts[, keep, drop = FALSE]
  ind  <- factor(meta$individual[keep])

  # pseudobulk: sum each gene across the cells of every individual
  # genes x individuals = counts %*% (cells x individuals 0/1 indicator)
  im <- sparse.model.matrix(~ 0 + ind)      # sparse-friendly aggregation
  pb <- as.matrix(cnt %*% im)
  colnames(pb) <- levels(ind)

  # one condition label per individual (columns of pb)
  group <- meta$condition[keep][match(levels(ind), ind)]

  # need exactly two conditions and enough individuals
  if (length(unique(group)) == 2 && ncol(pb) >= 4)
    results[[ct]] <- merit_de(pb, group)
}

# significant genes for one cell type, e.g.:
# subset(results[["CD4 T"]], padj < 0.05)
```

> If you use Seurat, `AggregateExpression(obj, group.by = c("individual","cell_type"))`
> produces the same per-individual pseudobulk counts; then run `merit_de()` per cell type
> exactly as above (pass the raw summed counts, not normalized data).

> Rule of thumb: at least 2–3 individuals per condition; more is better.
> 
### 5. Interpreting results

Table sorted by `padj`:

| Column | Meaning | Use |
|:--|:--|:--|
| `gene` | gene id | — |
| `stat` | calibrated z | sign = direction (positive = up in group 2); larger &#124;z&#124; = stronger |
| `pvalue` | raw p-value | usually not used directly |
| `padj` | BH-adjusted p (FDR) | **call significant at `padj < 0.05`** |

**FDR** = expected fraction of false calls among genes you declare significant.
MERIT-DE controls FDR strictly, giving a cleaner, more trustworthy significant-gene
list, while still detecting more true genes than Wilcoxon at small sample sizes.

---
