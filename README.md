# Comparative unproductive splicing across seven vertebrate species

Code and data for the comparative-splicing analysis: junction productivity classification,
differential splicing across developmental and adult-tissue contrasts, conservation of
unproductive cassette exons, and NMD-factor expression.

**[▶ Rendered figure notebook](https://bfairkun.github.io/20260825_comparativesplicing_paper/ComparativeSplicingFigures.html)**
— every panel, with the code that makes it. No setup required; the same file is committed at
[`docs/ComparativeSplicingFigures.html`](docs/ComparativeSplicingFigures.html).

## Reproducing the figures

The figure notebook reads only the small plot-ready tables in `output/` and the metadata in
`data/`, so it runs in seconds with no pipeline output present:

```bash
conda run -n py_general quarto render analysis/ComparativeSplicingFigures.qmd
```

It writes the individual data panels to `analysis/figures/`. The published figures are
multi-panel composites assembled by hand from those panels.

## Reproducing the data

```bash
cd code

# regenerate everything shipped in output/ from pipeline intermediates (cheap)
snakemake all_plot_data all_supplement_tables

# the whole pipeline, from downloading 1,993 FASTQs (weeks of compute, TBs of scratch)
snakemake all
```

Two targets are deliberately outside `all` because of their cost:

| Target | What it does |
|---|---|
| `all_dnds` | the 100-way multiz vertebrate alignment chain (~1.5 TB scratch) |
| `all_refresh_annotations` | re-fetch the Ensembl annotation tables via biomaRt |

`snakemake -n all` resolves the full DAG from nothing (~14.7k jobs) and is the check to run
after any edit.

## Layout

| Path | Contents |
|---|---|
| `analysis/` | the figure notebook and the panels it writes |
| `code/rules/` | project rules, including `PlotData.smk` (the terminal reduction rules) |
| `code/scripts/` | analysis and export scripts |
| `code/module_workflows/rna_seq/` | RNA-seq workflow, vendored as plain files |
| `code/config/` | sample, genome and contrast configuration |
| `data/` | shipped inputs: supplementary tables and pre-computed conservation data |
| `output/` | plot-ready tables read by the notebook |
| `code/scratch/` | disposable: parse-time derivatives and rule intermediates |

---

Secondary details — which differential-splicing implementation produced the distributed
results, vendored-code provenance, and known limitations — are in
[FOOTNOTES.md](FOOTNOTES.md).
