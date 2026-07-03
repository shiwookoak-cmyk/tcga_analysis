# TCGA Analysis Pipeline

Scripts for downloading TCGA data, running gene-wise survival analysis,
performing GSEA, visualizing KEGG metabolic pathways, calculating ssGSEA immune
signature scores, comparing mRNA expression, and creating oncoprints.

Downloaded TCGA data and generated outputs are saved under `data/` and are not
intended to be uploaded to GitHub.

## Repository Structure

```text
tcga_analysis/
|-- src/
|   |-- 01_download_tcga.R
|   |-- 02_gene_survival_analysis.R
|   |-- 03_run_gsea.R
|   |-- 04_map_kegg_pathways.R
|   |-- 05_run_ssgsea.R
|   |-- 06_analyze_expression.py
|   `-- 07_oncoprint.py
|-- renv/
|   |-- activate.R
|   `-- settings.json
|-- .Rprofile
|-- run_luad_pipeline.sh
|-- immune_signature.gmt
|-- renv.lock
|-- requirements.txt
`-- README.md
```

## Requirements

### R Packages

Exact R package versions are recorded in `renv.lock`. To restore the R
environment:

```r
install.packages("renv")
renv::restore()
```

The lockfile was generated with:

```text
R 4.6.0
Bioconductor 3.23
```

### Python Packages

```bash
pip install -r requirements.txt
```

## Example Workflow

Run the LUAD example pipeline:

```bash
bash run_luad_pipeline.sh
```

The script checks whether the expected LUAD input files already exist. If they
are present, the download step is skipped. Otherwise, TCGA data are downloaded
with `TCGAbiolinks`.

The run is logged to:

```text
data/logs/LUAD_pipeline_<timestamp>.log
```

## Individual Scripts

### 1. Download TCGA Data

```bash
Rscript src/01_download_tcga.R LUAD
```

Outputs:

```text
data/TCGA_LUAD/LUAD_rnaseq_counts.tsv
data/TCGA_LUAD/LUAD_rnaseq_log2tpm.tsv
data/TCGA_LUAD/LUAD_clinical.tsv
data/TCGA_LUAD/LUAD_snv.tsv
data/TCGA_LUAD/LUAD_cnv.tsv
```

### 2. Gene-Wise Survival Analysis

```bash
Rscript src/02_gene_survival_analysis.R LUAD median
Rscript src/02_gene_survival_analysis.R LUAD continuous
```

Outputs:

```text
data/LUAD_coxph_median.tsv
data/LUAD_coxph_continuous.tsv
```

### 3. GSEA

```bash
Rscript src/03_run_gsea.R data/LUAD_coxph_median.tsv
```

Outputs:

```text
data/LUAD_coxph_median_gsea_kegg.tsv
data/LUAD_coxph_median_gsea_reactome.tsv
data/LUAD_coxph_median_gsea_hallmark.tsv
data/LUAD_coxph_median_gsea_gobp.tsv
data/LUAD_coxph_median_gsea_gomf.tsv
data/LUAD_coxph_median_gsea_gocc.tsv
```

### 4. KEGG Metabolic Pathway Map

```bash
Rscript src/04_map_kegg_pathways.R data/LUAD_coxph_median_gsea_kegg.tsv
```

Output:

```text
data/LUAD_coxph_median_gsea_kegg_map.svg
```

### 5. ssGSEA Immune Signature Scoring

```bash
Rscript src/05_run_ssgsea.R LUAD immune_signature.gmt
```

Outputs:

```text
data/LUAD_ssgsea_scores.tsv
data/LUAD_ssgsea_signature_correlation_heatmap.svg
data/LUAD_ssgsea_sample_score_heatmap.svg
```

### 6. mRNA Expression Analysis

Edit the user settings in `main()` or import the `AnalyzeExpression` class into
a notebook:

```bash
python src/06_analyze_expression.py
```

### 7. Oncoprint

```bash
python src/07_oncoprint.py LUAD TP53 KRAS STK11 KEAP1 EGFR NF1
```

Outputs:

```text
data/LUAD_oncoprint_fisher.tsv
data/LUAD_oncoprint_fisher_heatmap.svg
data/LUAD_oncoprint.svg
```

## Notes

- TCGA data, GDC downloads, generated figures, intermediate R objects, and logs
  are excluded by `.gitignore`.
- `immune_signature.gmt` is used as an example GMT input for ssGSEA. Before
  publishing publicly, confirm that the file can be redistributed.
- Run scripts from the repository root so relative paths resolve correctly.
