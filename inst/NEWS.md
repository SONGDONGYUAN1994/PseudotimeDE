# Changes in PseudotimeDE

## Version 0.99.0 (First Bioconductor Submission)

### New Features
- Introduced `pseudotimeDE()`, a function for performing differential expression analysis along pseudotime using GAM models.
- Added `plotCurve()` for visualizing gene expression trends along pseudotime.
- Implemented `plotUncertainty()` to display pseudotime variability across subsamples.
- Included `runPseudotimeDE()`, a wrapper for running `pseudotimeDE()` on multiple genes efficiently.

### Significant User-Visible Changes
- Ensured compatibility with **SingleCellExperiment** and **Seurat** objects.
- Optimized computational performance with **parallel processing** (`BiocParallel`).




