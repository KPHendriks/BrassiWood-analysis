#!/usr/bin/env python3
"""
Summarize SH-aLRT, UFBoot, and branch lengths from IQ-TREE .treefile outputs.

New in this version:
- 'Passing all criteria' (count + %) added to overall_stats.txt
- Scatterplots use tiny points + compact marginal histograms (top/right ~1/5 size)
- Trendlines (NumPy) on ALL scatterplots, drawn on top with equation + R^2
- Branch-length scatter plots use log10(branch length) for x-axis
- All histograms (aLRT, UFBoot, branch length) include cumulative % line (right axis)

Typical usage:
  --input-glob "$RFIN/$project_name/3_results_iqtree_gene_trees/*/*_iqtree.treefile"
  --outdir     "$RFIN/$project_name/0_stats"
"""

import argparse, csv, glob, math, os, re, sys

try:
    from Bio import Phylo
except Exception:
    sys.stderr.write("[FATAL] Biopython is required (pip install biopython)\n")
    raise

HAS_PANDAS = False
HAS_MPL = False
HAS_NUMPY = False
try:
    import pandas as pd
    HAS_PANDAS = True
except Exception:
    pass

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except Exception:
    pass

try:
    import numpy as np
    HAS_NUMPY = True
except Exception:
    pass


def parse_args():
    ap = argparse.ArgumentParser(
        description="Summarize SH-aLRT, UFBoot, and branch lengths from IQ-TREE .treefile files."
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--input-glob", help='Glob, e.g. ".../3_results_iqtree_gene_trees/*/*_iqtree.treefile"')
    src.add_argument("--input-dir", help="Directory containing .treefile(s)")
    ap.add_argument("--pattern", default=".treefile", help='Suffix/pattern for files in --input-dir (default: ".treefile")')
    ap.add_argument("--outdir", required=True, help="Output directory (e.g., your 0_stats folder)")
    ap.add_argument("--alrt-threshold", type=float, default=80.0, help="Flag internal branches with SH-aLRT < this")
    ap.add_argument("--ufboot-threshold", type=float, default=95.0, help="Flag internal branches with UFBoot < this")
    ap.add_argument("--bl-threshold", type=float, default=1e-6, help="Flag branches with length < this")
    ap.add_argument("--gene-name-regex", default=r"^(.*)_iqtree\.treefile$", help="Regex to extract gene name from filename")
    ap.add_argument("--log10-bl-hist", action="store_true", help="Histogram of internal branch lengths on log10 scale")
    ap.add_argument("--pair-order", default="alrt,ufboot", choices=["alrt,ufboot", "ufboot,alrt"],
                    help="Order of two numbers in node labels when both present (default: alrt,ufboot)")
    # Kept for backward compatibility (ignored in this version):
    ap.add_argument("--scatter-style", default="auto", choices=["auto", "points", "hex"],
                    help="(Ignored) Previously controlled scatter vs hexbin")
    ap.add_argument("--hex-threshold", type=int, default=20000, help="(Ignored)")
    ap.add_argument("--hex-gridsize", type=int, default=60, help="(Ignored)")
    return ap.parse_args()


def find_files(args):
    if args.input_glob:
        return sorted(glob.glob(args.input_glob))
    else:
        return sorted(
            f for f in (os.path.join(args.input_dir, x) for x in os.listdir(args.input_dir))
            if (os.path.isfile(f) and args.pattern in os.path.basename(f))
        )


def extract_gene_name(basename, regex):
    m = re.match(regex, basename)
    if m:
        return m.group(1)
    return re.sub(r"\.treefile$", "", basename)


def extract_genome_from_path(path):
    # Expects .../3_results_iqtree_gene_trees/<GENOME>/<file>
    return os.path.basename(os.path.dirname(path))


PAIR_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)\s*$")
NUM_RE  = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*$")


def parse_supports(clade, pair_order):
    """
    Return (alrt, ufboot) floats or None.
    If only one number present, treat it as SH-aLRT (backward compatible).
    """
    alrt = None
    ufbt = None
    if clade.confidence is not None:
        try:
            alrt = float(clade.confidence)
            return alrt, ufbt
        except Exception:
            pass
    label = getattr(clade, "name", None)
    if not label:
        return alrt, ufbt
    m2 = PAIR_RE.match(label)
    if m2:
        a = float(m2.group(1)); b = float(m2.group(2))
        return (a, b) if pair_order == "alrt,ufboot" else (b, a)
    m1 = NUM_RE.match(label)
    if m1:
        alrt = float(m1.group(1))
        return alrt, ufbt
    return alrt, ufbt


def summarize_tree(tree_path, gene_name, genome, alrt_thr, ufbt_thr, bl_thr, pair_order):
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to parse {tree_path}: {e}\n")
        return [], None

    for i, clade in enumerate(tree.find_clades(order="preorder")):
        clade._node_id = i

    per_branch = []
    n_tips = 0
    alrts = []
    ufboots = []
    blens_internal = []
    n_internal = 0
    n_zero_alrt = 0
    n_low_alrt = 0
    n_low_ufbt = 0
    n_short_bl_internal = 0

    root = tree.root
    for clade in tree.find_clades(order="preorder"):
        is_leaf = clade.is_terminal()
        is_root = (clade is root)
        alrt, ufbt = (None, None)
        if not is_leaf:
            alrt, ufbt = parse_supports(clade, pair_order)

        blen = float(clade.branch_length) if (clade.branch_length is not None) else None

        if is_leaf:
            n_tips += 1
        else:
            n_internal += 1
            if alrt is not None:
                alrts.append(alrt)
                if alrt == 0:
                    n_zero_alrt += 1
                if alrt < alrt_thr:
                    n_low_alrt += 1
            if ufbt is not None:
                ufboots.append(ufbt)
                if ufbt < ufbt_thr:
                    n_low_ufbt += 1
            if blen is not None:
                blens_internal.append(blen)
                if blen < bl_thr:
                    n_short_bl_internal += 1

        per_branch.append(dict(
            genome=genome,
            gene=gene_name,
            node_id=clade._node_id,
            is_internal=(not is_leaf),
            is_root=is_root,
            support_alrt=(alrt if alrt is not None else ""),
            support_ufboot=(ufbt if ufbt is not None else ""),
            branch_length=(blen if blen is not None else ""),
            flag_low_alrt=int((not is_leaf) and (alrt is not None) and (alrt < alrt_thr)),
            flag_low_ufboot=int((not is_leaf) and (ufbt is not None) and (ufbt < ufbt_thr)),
            flag_short_bl=int(blen is not None and blen < bl_thr),
        ))

    def safe_mean(v): return sum(v) / len(v) if v else float("nan")
    def safe_median(v):
        if not v: return float("nan")
        s = sorted(v); k = len(s); m = k // 2
        return s[m] if k % 2 else 0.5 * (s[m-1] + s[m])

    per_gene_summary = dict(
        genome=genome,
        gene=gene_name,
        n_tips=n_tips,
        n_internal_branches=n_internal,
        n_internal_with_alrt=len(alrts),
        n_internal_with_ufboot=len(ufboots),
        n_zero_alrt=n_zero_alrt,
        frac_zero_alrt=(n_zero_alrt / n_internal if n_internal else float("nan")),
        n_low_alrt_lt_thr=n_low_alrt,
        frac_low_alrt_lt_thr=(n_low_alrt / n_internal if n_internal else float("nan")),
        n_low_ufboot_lt_thr=n_low_ufbt,
        frac_low_ufboot_lt_thr=(n_low_ufbt / n_internal if n_internal else float("nan")),
        n_short_bl_internal=n_short_bl_internal,
        frac_short_bl_internal=(n_short_bl_internal / n_internal if n_internal else float("nan")),
        mean_alrt_internal=safe_mean(alrts),
        median_alrt_internal=safe_median(alrts),
        mean_ufboot_internal=safe_mean(ufboots),
        median_ufboot_internal=safe_median(ufboots),
        mean_bl_internal=safe_mean(blens_internal),
        median_bl_internal=safe_median(blens_internal),
    )
    return per_branch, per_gene_summary


def write_csv(rows, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True
    )
    if HAS_PANDAS:
        pd.DataFrame(rows).to_csv(out_path, index=False)
    else:
        header = list(rows[0].keys()) if rows else []
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header); w.writeheader()
            for r in rows: w.writerow(r)


def _add_cumulative_line(ax_hist, counts, bin_edges):
    """Add cumulative % line on a secondary axis."""
    total = float(sum(counts)) if counts is not None else 0.0
    if total <= 0:
        return
    cum = []
    running = 0.0
    for c in counts:
        running += c
        cum.append(100.0 * running / total)
    ax2 = ax_hist.twinx()
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    ax2.plot(centers, cum, linewidth=1.2)
    ax2.set_ylabel("Cumulative (%)")
    ax2.set_ylim(0, 100)


def _scatter_with_marginals(x, y, xlabel, ylabel, title, vline=None, hline=None,
                            add_regression=True, out_png=None):
    """
    Render a main scatter with tiny points and marginal histograms (top/right).
    The marginal axes each take ~1/5 of the width/height.
    """
    if not HAS_MPL:
        return
    fig = plt.figure(figsize=(6.5, 6.0))
    gs = GridSpec(5, 5, figure=fig, height_ratios=[1, 0, 0, 0, 4], width_ratios=[4, 0, 0, 0, 1])
    # Main (bottom-left 4x4 block)
    ax_main = fig.add_subplot(gs[4, 0])
    # Top marginal (above main)
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    # Right marginal (to the right of main)
    ax_right = fig.add_subplot(gs[4, 4], sharey=ax_main)

    # Main scatter
    ax_main.scatter(x, y, s=2, alpha=0.25, zorder=1)
    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)
    ax_main.set_title(title)
    if vline is not None:
        ax_main.axvline(x=vline, linestyle="--", linewidth=1.2, zorder=2)
    if hline is not None:
        ax_main.axhline(y=hline, linestyle="--", linewidth=1.2, zorder=2)

    # Regression line (on top)
    if HAS_NUMPY and add_regression and len(x) >= 2:
        xx = np.asarray(x, dtype=float)
        yy = np.asarray(y, dtype=float)
        slope, intercept = np.polyfit(xx, yy, 1)
        xline = np.linspace(xx.min(), xx.max(), 200)
        yline = slope * xline + intercept
        ax_main.plot(xline, yline, linewidth=1.5, zorder=3)  # on top
        # R^2
        yhat = slope * xx + intercept
        ss_res = np.sum((yy - yhat) ** 2)
        ss_tot = np.sum((yy - yy.mean()) ** 2)
        r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else float("nan")
        ax_main.text(0.03, 0.97, f"y = {slope:.3f}x + {intercept:.3f}\nR² = {r2:.3f}",
                     transform=ax_main.transAxes, va="top")

    # Marginals: simple histograms
    ax_top.hist(x, bins=40)
    ax_right.hist(y, bins=40, orientation="horizontal")

    # Tidy marginal axes
    ax_top.tick_params(labelbottom=False)
    ax_right.tick_params(labelleft=False)
    for spine in ["top", "right"]:
        ax_main.spines[spine].set_visible(True)

    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=150)
    plt.close(fig)


def make_plots(branch_rows, outdir, log10_blen_hist, alrt_thr, ufbt_thr, bl_thr):
    if not HAS_MPL:
        sys.stderr.write("[INFO] matplotlib not available — skipping plots.\n")
        return

    # Internal branches with values present
    alrts  = [float(r["support_alrt"])  for r in branch_rows if r["is_internal"] and r["support_alrt"]  != ""]
    ufbt   = [float(r["support_ufboot"]) for r in branch_rows if r["is_internal"] and r["support_ufboot"] != ""]
    blens  = [float(r["branch_length"]) for r in branch_rows if r["is_internal"] and r["branch_length"] != ""]

    # SH-aLRT histogram with cumulative %
    if alrts:
        fig, ax = plt.subplots()
        counts, bins, _ = ax.hist(alrts, bins=40)
        ax.set_xlabel("SH-aLRT support (%) [internal branches]")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of SH-aLRT (internal branches)")
        if alrt_thr is not None:
            ax.axvline(x=alrt_thr, linestyle="--", linewidth=1.2)
        _add_cumulative_line(ax, counts, bins)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "overall_support_alrt_hist.png"), dpi=150)
        plt.close(fig)

    # UFBoot histogram with cumulative %
    if ufbt:
        fig, ax = plt.subplots()
        counts, bins, _ = ax.hist(ufbt, bins=40)
        ax.set_xlabel("UFBoot support (%) [internal branches]")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of UFBoot (internal branches)")
        if ufbt_thr is not None:
            ax.axvline(x=ufbt_thr, linestyle="--", linewidth=1.2)
        _add_cumulative_line(ax, counts, bins)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "overall_support_ufboot_hist.png"), dpi=150)
        plt.close(fig)

    # Branch length histogram with cumulative %
    if blens:
        data = [x for x in blens if x > 0] if log10_blen_hist else blens
        xlabel = "log10(Branch length) [internal]" if log10_blen_hist else "Branch length [internal]"
        if log10_blen_hist:
            data = [math.log10(x) for x in data]
        fig, ax = plt.subplots()
        counts, bins, _ = ax.hist(data, bins=40)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Count")
        ax.set_title("Distribution of branch lengths (internal)")
        thr = None
        if bl_thr is not None:
            thr = math.log10(bl_thr) if (log10_blen_hist and bl_thr > 0) else bl_thr
            ax.axvline(x=thr, linestyle="--", linewidth=1.2)
        _add_cumulative_line(ax, counts, bins)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "overall_branchlen_hist.png"), dpi=150)
        plt.close(fig)

    # Scatter: aLRT vs log10(branch length)
    if alrts and blens:
        xy = [(float(r["branch_length"]), float(r["support_alrt"])) for r in branch_rows
              if r["is_internal"] and r["branch_length"] != "" and r["support_alrt"] != "" and float(r["branch_length"]) > 0]
        if xy:
            x_raw, y = zip(*xy)
            x = [math.log10(v) for v in x_raw]
            _scatter_with_marginals(
                x, y,
                xlabel="log10(Branch length)",
                ylabel="SH-aLRT (%)",
                title="SH-aLRT vs log10(Branch length) [internal]",
                vline=(math.log10(bl_thr) if (bl_thr and bl_thr > 0) else None),
                hline=alrt_thr,
                add_regression=True,
                out_png=os.path.join(outdir, "scatter_support_alrt_vs_branchlen.png")
            )

    # Scatter: UFBoot vs log10(branch length)
    if ufbt and blens:
        xy = [(float(r["branch_length"]), float(r["support_ufboot"])) for r in branch_rows
              if r["is_internal"] and r["branch_length"] != "" and r["support_ufboot"] != "" and float(r["branch_length"]) > 0]
        if xy:
            x_raw, y = zip(*xy)
            x = [math.log10(v) for v in x_raw]
            _scatter_with_marginals(
                x, y,
                xlabel="log10(Branch length)",
                ylabel="UFBoot (%)",
                title="UFBoot vs log10(Branch length) [internal]",
                vline=(math.log10(bl_thr) if (bl_thr and bl_thr > 0) else None),
                hline=ufbt_thr,
                add_regression=True,
                out_png=os.path.join(outdir, "scatter_support_ufboot_vs_branchlen.png")
            )

    # Scatter: SH-aLRT vs UFBoot (with regression)
    if alrts and ufbt:
        xy = [(float(r["support_alrt"]), float(r["support_ufboot"])) for r in branch_rows
              if r["is_internal"] and r["support_alrt"] != "" and r["support_ufboot"] != ""]
        if xy:
            x, y = zip(*xy)
            _scatter_with_marginals(
                list(x), list(y),
                xlabel="SH-aLRT (%)",
                ylabel="UFBoot (%)",
                title="SH-aLRT vs UFBoot (internal branches)",
                vline=alrt_thr,
                hline=ufbt_thr,
                add_regression=True,
                out_png=os.path.join(outdir, "scatter_alrt_vs_ufboot.png")
            )


def main():
    args = parse_args()
    files = find_files(args)
    if not files:
        sys.stderr.write("[ERROR] No .treefile inputs found.\n")
        sys.exit(2)

    os.makedirs(args.outdir, exist_ok=True)

    branch_rows, per_gene_rows = [], []
    for fp in files:
        bn = os.path.basename(fp)
        gene = extract_gene_name(bn, args.gene_name_regex)
        genome = extract_genome_from_path(fp)
        per_branch, per_gene = summarize_tree(
            fp, gene, genome,
            args.alrt_threshold, args.ufboot_threshold, args.bl_threshold,
            args.pair_order
        )
        if per_branch and per_gene:
            branch_rows.extend(per_branch)
            per_gene_rows.append(per_gene)

    if not per_gene_rows:
        sys.stderr.write("[ERROR] No trees could be summarized.\n")
        sys.exit(3)

    write_csv(branch_rows, os.path.join(args.outdir, "branch_level.csv"))
    write_csv(per_gene_rows, os.path.join(args.outdir, "per_gene_summary.csv"))

    # Totals
    total_internal = sum(r["n_internal_branches"] for r in per_gene_rows if r["n_internal_branches"] == r["n_internal_branches"])
    total_zero = sum(r["n_zero_alrt"] for r in per_gene_rows if r["n_zero_alrt"] == r["n_zero_alrt"])
    total_low_alrt = sum(r["n_low_alrt_lt_thr"] for r in per_gene_rows if r["n_low_alrt_lt_thr"] == r["n_low_alrt_lt_thr"])
    total_low_ufbt = sum(r["n_low_ufboot_lt_thr"] for r in per_gene_rows if r["n_low_ufboot_lt_thr"] == r["n_low_ufboot_lt_thr"])
    total_short = sum(r["n_short_bl_internal"] for r in per_gene_rows if r["n_short_bl_internal"] == r["n_short_bl_internal"])

    # Passing ALL criteria (conservative: missing values fail)
    passing_all = 0
    for r in branch_rows:
        if not r["is_internal"]:
            continue
        try:
            a = float(r["support_alrt"]) if r["support_alrt"] != "" else None
            b = float(r["support_ufboot"]) if r["support_ufboot"] != "" else None
            bl = float(r["branch_length"]) if r["branch_length"] != "" else None
        except Exception:
            a = b = bl = None
        ok = (
            a is not None and a >= args.alrt_threshold and
            b is not None and b >= args.ufboot_threshold and
            bl is not None and bl >= args.bl_threshold
        )
        if ok:
            passing_all += 1

    with open(os.path.join(args.outdir, "overall_stats.txt"), "w") as fh:
        fh.write(f"Files parsed: {len(per_gene_rows)}\n")
        fh.write(f"Internal branches (sum): {total_internal}\n")
        pct_zero = (total_zero / total_internal * 100.0) if total_internal else float('nan')
        fh.write(f"Zero aLRT (sum): {total_zero} ({pct_zero:.2f}%)\n")
        pct_low_alrt = (total_low_alrt / total_internal * 100.0) if total_internal else float('nan')
        fh.write(f"Low SH-aLRT < {args.alrt_threshold} (sum): {total_low_alrt} ({pct_low_alrt:.2f}%)\n")
        pct_low_ufbt = (total_low_ufbt / total_internal * 100.0) if total_internal else float('nan')
        fh.write(f"Low UFBoot < {args.ufboot_threshold} (sum): {total_low_ufbt} ({pct_low_ufbt:.2f}%)\n")
        pct_short = (total_short / total_internal * 100.0) if total_internal else float('nan')
        fh.write(f"Short internal branches < {args.bl_threshold} (sum): {total_short} ({pct_short:.2f}%)\n")
        pct_pass = (passing_all / total_internal * 100.0) if total_internal else float('nan')
        fh.write(f"Passing all criteria (≥aLRT {args.alrt_threshold}, ≥UFBoot {args.ufboot_threshold}, ≥BL {args.bl_threshold}): "
                 f"{passing_all} ({pct_pass:.2f}%)\n")

    make_plots(
        branch_rows,
        args.outdir,
        log10_blen_hist=args.log10_bl_hist,
        alrt_thr=args.alrt_threshold,
        ufbt_thr=args.ufboot_threshold,
        bl_thr=args.bl_threshold
    )

    print("[OK] Wrote into:", args.outdir)
    print(" - branch_level.csv")
    print(" - per_gene_summary.csv")
    print(" - overall_stats.txt")
    if HAS_MPL:
        print(" - overall_support_alrt_hist.png")
        print(" - overall_support_ufboot_hist.png")
        print(" - overall_branchlen_hist.png")
        print(" - scatter_support_alrt_vs_branchlen.png")
        print(" - scatter_support_ufboot_vs_branchlen.png")
        print(" - scatter_alrt_vs_ufboot.png")


if __name__ == "__main__":
    main()
