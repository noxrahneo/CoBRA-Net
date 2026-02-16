# CoBRA-Net (Thesis Workspace)

This repository is the organized workspace for my thesis final version. It will evolve as the project progresses.

**Last updated**: February 16, 2026

## Current Structure

- `scripts/`: Main analysis workflow scripts (.py). These cover environment setup, data loading, preprocessing, QC analysis, and network building.
- `utils/`: Helper utilities used by the main scripts.

## Workflow Outline

1. Set up environment and load scRNA-seq data.
2. Preprocess data (basic scRNA-seq pipeline).
3. Initial data quality assessment (pre-QC) and visualization.
4. Perform QC and filtering (retain required data only).
5. Repeat quality assessment (post-QC) to show improvements and generate plots for the thesis.
6. Normalize data and select highly variable genes.
7. Annotate cell types.
8. Merge by cancer type and extract epithelial-only subsets:
	 - 1 healthy, 1 preneoplastic, 4 cancerous.
9. Compute correlations.
10. Build networks (individual and combined) to study cancer development.
11. Identify hub genes.
12. Downstream analyses: enrichment analysis, expression directionality,
		module comparison and therapeutic targeting to identify candidate drugs.

## Notes

- QC/preprocessing follows standard scRNA-seq best practices:
	https://www.sc-best-practices.org/introduction/analysis_tools.html
- Future additions may include garbage-collection helpers and expanded tooling.

## Automation

This repo includes a `hooks/` folder to share git hook scripts. To enable them locally:

```bash
git config core.hooksPath hooks
```

Current hooks:

- `pre-commit`: updates the **Last updated** date in README.md.
