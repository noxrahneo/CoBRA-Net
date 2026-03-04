# Interactive network filtering: biological and computational meaning

This document explains what the sliders in the interactive network HTML do, and how to interpret them for biology.

## What the graph already contains before filtering

The interactive graph is generated from the **precomputed sparse network** for a condition (e.g., `Normal`).

- Nodes = genes kept by upstream network construction.
- Edges = gene-gene connections kept by upstream sparsification and edge selection.
- Node metadata includes weighted degree (`wd`), i.e., summed edge strengths connected to that gene.
- Edge metadata includes edge weight (`ew`), i.e., connection strength for that pair.

Important: the slider filtering does **not** recompute correlations or network inference. It only subsets what is already in this precomputed graph.

## What “Genes shown” does

The **Genes shown** slider ranks all currently available nodes by weighted degree (`wd`) and keeps the top percentage.

- Higher %: includes more genes, including lower-connectivity genes.
- Lower %: focuses on hub-like genes with strongest overall connectivity in this graph.

Biological interpretation:

- Decreasing this slider enriches for genes with stronger network centrality (putative coordinators/hubs in the inferred module structure).
- Increasing it reveals more peripheral genes that may be condition-relevant but less connected.

## What “Connections shown” does

After gene filtering, remaining candidate edges are ranked by edge weight (`ew`) and only the top percentage is kept.

- Higher %: denser network, more medium/weak retained links.
- Lower %: emphasizes strongest co-expression links among the visible genes.

Biological interpretation:

- Decreasing this slider highlights the strongest pairwise associations.
- Increasing it restores broader local context around those associations.

## What “adding/removing genes or connections” means here

In this interactive view, adding/removing is **visual inclusion/exclusion only**.

- It does not create or destroy biological relationships.
- It does not rerun the solver or change node coordinates.
- It does not alter the saved upstream results/tables.

So, slider changes are best interpreted as **different visibility thresholds** over the same inferred network.

## Why layout stays stable

The network physics is disabled (`enabled: false`) for this view. Filtering updates the visible node/edge sets but does not trigger a new force-layout solve.

This allows direct visual comparison across slider settings without layout drift.

## Export behavior

Two high-resolution PNG export modes are provided:

- `Save view (PNG, high-res)`: exports current canvas view as-is.
- `Save view (PNG white bg, high-res)`: exports same view with an opaque white background.

Both export filenames include:

- visible gene count
- visible connection count
- detected solver name

This helps keep exported figures traceable to the exact view state.
