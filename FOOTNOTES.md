# Footnotes

Details that would clutter the [README](README.md): which differential-splicing
implementation produced the distributed results, where the vendored code came from, and what
this bundle does not reproduce. Nothing here is needed to render the figures or run the
workflow.

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

## Junction classification: which leafcutter2 produced what

Upstream leafcutter2 commit `16b5603` (2026-04-10, "recursive splicing fixes") changed ORF
path-finding in `solve_NMD` so that a junction abutting the previous one can extend a path
(`last_pos == j0 + 1`), and so that zero-length exons are skipped during translation. It is
unconditional — there is no runtime flag — and it strictly widens the set of reachable ORF
paths, so classifications move: `Coding` flips only False→True, `UTR` only True→False, and
after `UTR & !Coding → NE` reconciliation the class transitions are overwhelmingly UP→PR and
NE→PR. No junction is gained or lost, so diffing junction *sets* shows nothing. (Not to be
confused with `--remove-junctions`, commit `2a2094c`, whose message also mentions recursive
splicing but which is optional, off by default, and a different mechanism.)

**The classifications distributed here predate that change**, on two independent lines of
evidence:

- The vendored classifier is the older script-based
  `code/scripts/leafcutter2_daiuc/scripts/ForwardSpliceJunctionClassifier.py`
  (daiuc/leafcutter2 `eb99084`), not the restructured `src/leafcutter2/classifier.py` the fix
  applies to. Its path chaining is the pre-fix form —
  `if strand == "+" and last_pos < j0 and j0-last_pos < exonLcutoff` — with no abutting-junction
  clause and no zero-length-exon skip. The second-stage "extra UTR rule"
  (`if utr and j[0] in coding_5_prime and ...`) is present and identical, which is the route by
  which a junction's `UTR` can flip even when its own `Coding` does not: `coding_5_prime` is
  built from junctions that passed `solve_NMD`, so a sibling sharing the 5' splice site is
  enough.
- The `_junction_classifications.txt` files that the shipped tallies were reduced from are dated
  **2024-10-24**, about eighteen months before the fix.

**This does not need a version knob here.** leafcutter2 is vendored as plain tracked files and
invoked by path (`python scripts/leafcutter2_daiuc/scripts/...`); nothing installs it from a
package or a moving ref. The workflow therefore reproduces its own classifications regardless of
what upstream does, which is the property a floating install would not have. The vendored copy
is md5-identical (`f6f7bfc3…`) to the copy in the working repository that actually produced the
published files.

What this does mean is that these classifications will **not** match output from current
leafcutter2 (v2.0.1, `7dedb0a`), and will differ in the direction described above. The
*magnitude* here is unmeasured — quantifying it would mean re-running classification for all
seven species, which is outside what this bundle re-runs. For scale, an independent simulation
benchmark measured `Coding` flips at 0.04–0.06% of junction-condition rows with every metric
median unmoved, but that was different annotations and different data, so treat it as an order
of magnitude and not a prediction.

One related hazard that does *not* apply here: leafcutter2 v2.0.1 fails on fresh installs
because bedparse 0.2.3 imports `pkg_resources`, removed in setuptools 81, surfacing as a
misleading `'NoneType' object has no attribute 'splitlines'`. Both bedparse envs in this
project already pin `setuptools==63.4.1`, so the failure cannot occur.

Credit: the upstream change, its observable signature, and the byte-for-byte bisection that
identified `2c161e2` as the last pre-fix commit were worked out in the
`20260825_Leaf2Simulation_PublicationFigs` project; see its `FOOTNOTES.md`, section
"leafcutter2: which version produced what".

## Juvenile versus adult testis: two group definitions

> **The supplementary *UPF3A*/*UPF3B* panel as published grouped testis libraries by the
> onset-of-meiosis developmental landmark, while Fig 3c grouped them by the supplementary
> contrast table. The figure notebook now uses the contrast table for both.**

There are two ways to decide whether a testis library is juvenile or adult. The landmark rule
calls a library adult once its within-species percent-rank stage reaches the onset of meiosis
in spermatogenesis. The contrast table (`data/20260820_SupplementalTable_Contrasts.tsv`) names
the two arms explicitly, and is what every differential-splicing contrast in the paper was run
on, Fig 3c included.

The two are not independent. The contrast groups were originally drawn *by* the landmark rule,
from the same ordinal stage table, and were consistent with it. They were then corrected: the
libraries sitting exactly on the landmark were inspected, and in some species those borderline
samples looked pre-meiotic while in others they looked post-meiotic, so the boundary stage was
assigned per species on that evidence. The contrast table is therefore the landmark rule plus a
hand adjudication of the one ambiguous stage.

The two groupings consequently partition the same libraries — per-species totals are equal in
all seven species — and disagree about exactly the libraries that adjudication moved: 11
libraries in five species, every one of them sitting precisely on the onset-of-meiosis
landmark, where the rule's `>=` resolves the tie toward adult for everybody. Rabbit and chicken
have a library at that stage too, but there the borderline sample was judged post-meiotic, so
the table agrees with the rule — which is why those are the two species whose group sizes
already matched.

The contrast table is the grouping to prefer: it carries information the stage table does not,
namely what the borderline libraries actually look like.

Consequences are confined to rat, which goes from non-significant under the landmark grouping
(P = 0.1, 0.3, 0.2 for *UPF3A*, *UPF3B*, ratio) to significant under the contrast table
(P = 4e-05, 0.02, 2e-04). No effect changes direction in any species.

**[▶ Rendered comparison notebook](https://bfairkun.github.io/20260825_comparativesplicing_paper/NMD_GroupingDiscrepancy.html)**
— draws the panel both ways from the shipped tables, plots every library against its ordinal
stage with the developmental landmarks marked, and names the libraries that move. Committed at
[`docs/NMD_GroupingDiscrepancy.html`](docs/NMD_GroupingDiscrepancy.html); rebuild with
`conda run -n py_general quarto render analysis/NMD_GroupingDiscrepancy.qmd`.

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
