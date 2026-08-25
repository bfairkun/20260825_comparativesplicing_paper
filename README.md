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
| `code/module_workflows/rna_seq/` | vendored RNA-seq workflow (see below) |
| `code/config/` | sample, genome and contrast configuration |
| `data/` | shipped inputs: supplementary tables and pre-computed conservation data |
| `output/` | plot-ready tables read by the notebook |
| `code/scratch/` | disposable: parse-time derivatives and rule intermediates |

## Differential splicing: which implementation produced what

> **The differential-splicing results distributed here were generated with the original
> `leafcutter_ds.R`; the Python `leafcutter-ds` implementation that the workflow now calls was
> not available at the time of manuscript preparation.**

The workflow calls the Python implementation
([leafcutter2/leafcutter-ds](https://github.com/leafcutter2/leafcutter-ds)) through
`code/scripts/leafcutter_ds_py.py`, a thin wrapper that reconciles two differences from the R
tool so nothing downstream had to change:

- **Baseline group.** `leafcutter-ds` defaults `--baseline_group` to the literal string
  `Control`. When no group has that name it still picks a baseline, which flips the sign of
  every `deltapsi` and emits a `psi_Control` column matching no input group. The wrapper
  passes the baseline explicitly, using the same rule as `leafcutter_ds.R`: the group label
  that appears **first in the groups file**. `leafcutter_ds.R` sets its factor levels with
  `unique(meta$group)`, i.e. order of first appearance, and sorts only if that column parses
  as numeric — so group-file row order is load-bearing and fixes the sign of `deltapsi`.
  The group files are generated in the order actually used (see below).
- **Table schema.** The Python tool suffixes columns by group (`logef_<g2>`, `psi_<g1>`,
  `deltapsi_<g2>`) and puts `cluster` first in the significance table with an extra `genes`
  column. The wrapper renames and reorders both tables to the R schema.

Agreement on a 400-cluster test against the archived R output for
`Neonate_vs_Adult.Liver.Chicken`: `psi` r=1.0000, `deltapsi` r=0.9976 (median |diff|
1.3e-4), `loglr` r=0.9979, `p` r=0.9963, `status` 100%, `df` exact, and every intron id
matched by exact string with no coordinate shift. `logef` is the one column that genuinely
diverges (r=0.89) because the two implementations regularise the fit differently; nothing
downstream reads it, and it is carried into the tidy tables for reference only.

Installation note: the PyPI wheel (`leafcutter==2.0.2`) is broken — it contains only
dist-info and the two console scripts, with no `leafcutter` package — so
`code/envs/leafcutter_ds_py.yml` installs from git, pinned to commit `bdf92d37`, the version
the agreement above was measured against. The pin is deliberate: this implementation is new
and liable to be patched, so moving it means re-running that comparison. That env pulls torch
and the CUDA runtime via `pyro-ppl` and is roughly 3 GB.

**Group files and their row order.** Because the baseline is the first label in the file,
`rule MakeContrastGroupFiles` must reproduce row order, not just the sample→group mapping. It
does: all 223 files it generates from `data/20260820_SupplementalTable_Contrasts.tsv` are
byte-identical, order included, to the ones used for the published analysis. For every one of
those contrasts the first-in-file and alphabetically-first labels happen to coincide — the
developmental contrasts are all `Early`/`Late`, and the adult-tissue pairs are enumerated from
an alphabetically ordered tissue list — so the two conventions cannot be distinguished from
the shipped results. They would diverge for a contrast added with a different row order.

**Clustering** is unchanged: the original leafcutter's `leafcutter_cluster_regtools.py`,
which is the correct pairing with the regtools junction extraction used upstream (it applies
the regtools block-offset correction). Note the coordinate convention that follows from it:
the leafcutter cluster key's `stop` is one past the intron end, so
`tidy_leafcutter_ds_results.R` forms `Intron_coord` as `{chrom}:{start}-{stop-1}` to line up
with the leafcutter2 classifications.

## Vendored code

No git submodules; the following are vendored as plain files at these commits.

| Path | Source | Commit |
|---|---|---|
| `code/module_workflows/rna_seq` | bfairkun/snakemake-workflow_rna-seq | `517f927` |
| `code/module_workflows/rna_seq/scripts/leafcutter` | davidaknowles/leafcutter | `2c9907e` |
| `code/scripts/leafcutter2_bfairkun` | bfairkun/leafcutter2 | `2bc249f` |
| `code/scripts/leafcutter2_daiuc` | daiuc/leafcutter2 | `eb99084` |

Installed rather than vendored, but pinned for the same reason:

| Env | Source | Commit |
|---|---|---|
| `code/envs/leafcutter_ds_py.yml` | leafcutter2/leafcutter-ds | `bdf92d37` |

The rna-seq workflow is pinned to the commit the parent repository recorded, not the newer
commit that happened to be checked out, which changes the sample-sheet schema and adds
per-sample network calls at parse time.

## Known limitations

- The dN/dS and PhyloP chain (`all_dnds`) is the one part of the workflow not exercised by
  the dry-run check, because its small end products are shipped in `data/conservation/`
  rather than rebuilt.
- FASTQ download uses aspera when `aspera_key` is set in `config/config.yaml` and plain FTP
  otherwise. Only the FTP path is portable.
