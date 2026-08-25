#!/usr/bin/env python3
"""Run the Python leafcutter-ds and emit output in the original leafcutter_ds.R schema.

The Python implementation (https://github.com/leafcutter2/leafcutter-ds) takes the same
inputs as leafcutter_ds.R and writes the same two filenames, but its table schemas differ,
so downstream code would break on a straight swap. This wrapper handles the two differences:

1. **Baseline group.** leafcutter-ds needs `-0/--baseline_group` naming one of the groups.
   Its default is the literal string "Control"; if no group has that name it silently picks
   a baseline anyway, which FLIPS THE SIGN of every deltapsi and emits a column called
   `psi_Control` that does not correspond to any group in the input. leafcutter_ds.R instead
   takes the first level of `factor(groups, unique(groups))`, i.e. the label appearing first
   in the groups file (see `baseline_group` below). This wrapper derives that label and
   passes it explicitly, so deltapsi keeps the R sign convention
   (deltapsi = psi[second group] - psi[baseline]).

2. **Column names and order.**

       effect_sizes           R:  intron  logef        <g1>      <g2>     deltapsi
                             Py:  intron  logef_<g2>   psi_<g1>  psi_<g2> deltapsi_<g2>

       cluster_significance   R:  status  loglr  df  p  cluster  p.adjust
                             Py:  cluster status loglr df p      p.adjust  genes

   The wrapper renames and reorders to the R schema and drops the extra `genes` column, so
   `tidy_leafcutter_ds_results.R` and everything downstream is unchanged.

Numerical agreement on a 400-cluster test against the archived R output for
Neonate_vs_Adult.Liver.Chicken: psi r=1.0000, deltapsi r=0.9976 (median |diff| 1.3e-4),
loglr r=0.9979, p r=0.9963, status 100%, df exact, and every intron id matched by exact
string (no coordinate shift). `logef` is the one column that genuinely diverges (r=0.89) --
the two implementations regularise the fit differently. Nothing downstream reads `logef`;
it is passed through into the tidy tables for reference only.

Tested against leafcutter-ds commit bdf92d377f5dfbe5d643b5c217dc7047d70742e8 (2026-03-30),
which is the commit pinned in ../envs/leafcutter_ds_py.yml. The Python implementation is new
and may change; if that pin is moved, re-run the comparison above before trusting the output.
"""
import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd

R_CLUSTER_COLS = ["status", "loglr", "df", "p", "cluster", "p.adjust"]


def baseline_group(groups_file):
    """The group leafcutter_ds.R treats as baseline: the first label in FILE ORDER.

    leafcutter_ds.R builds its factor levels as

        group_names = unique(meta$group)              # order of first appearance
        if (is.numeric(meta$group)) group_names = sort(group_names)
        meta$group = factor(meta$group, group_names)

    so the baseline is whichever label appears first in the groups file -- NOT the
    alphabetically first one. R only sorts when the group column parses as numeric, which
    never happens for the string labels this project uses. `pandas.unique` preserves order
    of appearance, matching R's `unique`; the row order of the groups file is therefore
    load-bearing, since it fixes the sign of every deltapsi.
    """
    g = pd.read_csv(groups_file, sep=r"\s+", header=None, usecols=[0, 1],
                    names=["sample", "group"], dtype=str)
    groups = list(pd.unique(g.group))
    if len(groups) != 2:
        sys.exit(f"expected exactly 2 groups in {groups_file}, found {groups}")
    return groups[0], groups[1]


def rewrite_effect_sizes(path, g1, g2):
    d = pd.read_csv(path, sep="\t")
    ren = {f"logef_{g2}": "logef", f"psi_{g1}": g1, f"psi_{g2}": g2,
           f"deltapsi_{g2}": "deltapsi"}
    missing = [k for k in ren if k not in d.columns]
    if missing:
        sys.exit(f"leafcutter-ds effect_sizes lacks expected columns {missing}; "
                 f"got {list(d.columns)}. Was --baseline_group set to {g1!r}?")
    d = d.rename(columns=ren)
    d[["intron", "logef", g1, g2, "deltapsi"]].to_csv(path, sep="\t", index=False)


def rewrite_cluster_significance(path):
    d = pd.read_csv(path, sep="\t")
    missing = [c for c in R_CLUSTER_COLS if c not in d.columns]
    if missing:
        sys.exit(f"leafcutter-ds cluster_significance lacks {missing}; got {list(d.columns)}")
    d[R_CLUSTER_COLS].to_csv(path, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("counts_file")
    ap.add_argument("groups_file")
    ap.add_argument("-o", "--output_prefix", default="leafcutter_ds")
    ap.add_argument("-i", "--min_samples_per_intron", default=None)
    ap.add_argument("-g", "--min_samples_per_group", default=None)
    ap.add_argument("-c", "--min_coverage", default=None)
    ap.add_argument("-s", "--max_cluster_size", default=None)
    ap.add_argument("-e", "--exon_file", default=None)
    ap.add_argument("-p", "--num_threads", default=None)
    args = ap.parse_args()

    exe = shutil.which("leafcutter-ds")
    if exe is None:
        sys.exit("leafcutter-ds not on PATH. NOTE: the PyPI wheel (leafcutter==2.0.2) is "
                 "broken -- it ships no package, only console scripts -- so install from "
                 "source: pip install git+https://github.com/leafcutter2/leafcutter-ds.git")

    g1, g2 = baseline_group(args.groups_file)
    cmd = [exe, "-o", args.output_prefix, "-0", g1]
    for flag, val in [("-i", args.min_samples_per_intron),
                      ("-g", args.min_samples_per_group),
                      ("-c", args.min_coverage),
                      ("-s", args.max_cluster_size),
                      ("-e", args.exon_file),
                      ("-p", args.num_threads)]:
        if val is not None:
            cmd += [flag, str(val)]
    cmd += [args.counts_file, args.groups_file]
    print(f"baseline group (R factor-level convention): {g1!r}; "
          f"deltapsi = psi[{g2}] - psi[{g1}]", file=sys.stderr)
    print("+ " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

    rewrite_effect_sizes(f"{args.output_prefix}_effect_sizes.txt", g1, g2)
    rewrite_cluster_significance(f"{args.output_prefix}_cluster_significance.txt")
    print("rewrote both tables into the leafcutter_ds.R schema", file=sys.stderr)


if __name__ == "__main__":
    main()
