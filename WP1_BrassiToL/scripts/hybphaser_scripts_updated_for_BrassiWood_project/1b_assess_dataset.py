#!/usr/bin/env python3
# 1b_assess_dataset.py  (v4 — adds per-cell min length fraction; 2025-09-24)
#
# New in v4:
#   - --min_len_frac_cell (float, default 0.0 = off). If > 0, for each (locus, sample),
#     we require consensus_length / target_length >= threshold. Otherwise that cell is
#     treated as missing (NA) in BOTH SNP and length tables before any other filtering.
#
# Still included from previous revisions:
#   - Correct n_loci_raw / nloci_remaining math (computed from pre-clean presence + union of flags).
#   - Scatter points without borders (edgecolors='none').
#   - "Failed sample" requires SNPs all NA AND lengths all NA or zero.
#   - HybPiper TXT scan cached to CSV when provided.
#
import argparse
from pathlib import Path
import sys
import re
from collections import defaultdict, OrderedDict
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def read_csv_or_die(p: Path, what: str) -> pd.DataFrame:
    if not p.exists():
        sys.exit(f"[FATAL] Missing {what}: {p}")
    return pd.read_csv(p, index_col=0, low_memory=False)


def infer_locus_from_text_line(s: str) -> str:
    s = re.sub(r"\s.*$", "", str(s).strip())
    return s.rsplit("-", 1)[-1] if "-" in s else s

def scan_hybpiper_txt(hdir: Path, samples: List[str], genomes: Optional[List[str]]) -> List[Tuple[str, str, str]]:
    out: List[Tuple[str, str, str]] = []
    if hdir is None or not hdir.exists():
        return out
    for s in samples:
        if genomes:
            candidates = []
            for g in genomes:
                fname = f"{s}_{g}_genes_with_long_paralog_warnings.txt"
                candidates.extend(hdir.rglob(fname))
        else:
            candidates = list(hdir.rglob(f"{s}_*_genes_with_long_paralog_warnings.txt"))
        for fp in candidates:
            try:
                for ln in fp.read_text().splitlines():
                    ln = ln.strip()
                    if not ln or ln.startswith("#"):
                        continue
                    m = re.match(rf"^{re.escape(s)}_([^_]+)_genes_with_long_paralog_warnings\.txt$", fp.name)
                    gname = m.group(1) if m else ""
                    out.append((s, gname, infer_locus_from_text_line(ln)))
            except Exception:
                continue
    return out


def parse_targets_arg(items):
    d = OrderedDict()
    for it in items:
        if "=" not in it:
            sys.exit(f"[FATAL] --targets entry must be GENOME=/path/to.fasta, got '{it}'")
        k, v = it.split("=", 1)
        k = k.strip(); v = v.strip()
        if not k or not v:
            sys.exit(f"[FATAL] malformed --targets entry: '{it}'")
        d[k] = Path(v)
    return d

def gene_key_from_header(header: str, mode: str) -> str:
    if mode in ("dash", "full"):
        rid = re.sub(r"\s.*$", "", header)
    elif mode == "space":
        rid = header.split()[0]
    else:
        rid = re.sub(r"\s.*$", "", header)
    if mode == "dash":
        return rid.rsplit("-", 1)[-1] if "-" in rid else rid
    elif mode == "full":
        return rid
    elif mode == "space":
        return rid.rsplit("-", 1)[-1] if "-" in rid else rid
    return rid

def parse_targets_lengths(target_map: dict, fmt: str, id_mode: str):
    is_aa = fmt.upper() == "AA"
    max_per_gene = {}
    counts = {}
    gene_to_genome = {}
    per_genome_keyset = {g:set() for g in target_map.keys()}

    for genome, fasta_path in target_map.items():
        if not fasta_path.exists():
            print(f"[WARN] Targets file missing for genome {genome}: {fasta_path}", file=sys.stderr)
            counts[genome] = (0, 0)
            continue
        seen = {}
        nrecs = 0
        try:
            for rec in SeqIO.parse(str(fasta_path), "fasta"):
                nrecs += 1
                key = gene_key_from_header(rec.id, id_mode)
                L = len(rec.seq) * (3 if is_aa else 1)
                if key not in seen or L > seen[key]:
                    seen[key] = L
                if key not in max_per_gene or L > max_per_gene[key]:
                    max_per_gene[key] = L
                if key not in gene_to_genome:
                    gene_to_genome[key] = genome
                per_genome_keyset[genome].add(key)
        except Exception as e:
            sys.exit(f"[FATAL] Failed reading targets {fasta_path}: {e}")
        counts[genome] = (nrecs, len(seen))
    return max_per_gene, counts, gene_to_genome, per_genome_keyset


def iqr_threshold(values: np.ndarray) -> float:
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return np.nan
    q1, q3 = np.quantile(v, [0.25, 0.75])
    return q3 + 1.5 * (q3 - q1)


SCATTER_ALPHA = 0.5
SCATTER_S_SMALL = 10
SCATTER_S_MED = 18

def maybe_pdf_png(base: Path, plot_fn, size=(12, 8), dpi=160, make_pdf=False):
    plt.figure(figsize=size); plot_fn(); plt.tight_layout()
    plt.savefig(base.with_suffix(".png"), dpi=dpi); plt.close()
    if make_pdf:
        plt.figure(figsize=size); plot_fn(); plt.tight_layout()
        plt.savefig(base.with_suffix(".pdf")); plt.close()

def color_cycle(n):
    from matplotlib.cm import get_cmap
    cmap = get_cmap("tab10")
    return [cmap(i % 10) for i in range(n)]


def order_loci_by_genome_then_name(loci_iter, gene_to_genome_map):
    def keyfunc(loc):
        return (str(gene_to_genome_map.get(loc, "")), str(loc))
    return sorted(map(str, loci_iter), key=keyfunc)


def main():
    ap = argparse.ArgumentParser(description="HybPhaser 1b: assess missing data, remove paralogs, compute LH & AD.")
    ap.add_argument("--out_root", required=True, type=Path)
    ap.add_argument("--subset_name", default="")
    ap.add_argument("--targets", action="append", default=[],
                    help="Repeatable: GENOME=/path/to/targets.fasta (merged by gene key with max length).")
    ap.add_argument("--targets_format", choices=["DNA", "AA"], default="DNA")
    ap.add_argument("--target_id_mode", choices=["dash", "space", "full"], default="dash")
    ap.add_argument("--thr_samples_prop_loci", type=float, default=0.2)
    ap.add_argument("--thr_samples_prop_target", type=float, default=0.2)
    ap.add_argument("--thr_loci_prop_samples", type=float, default=0.2)
    ap.add_argument("--thr_loci_prop_target", type=float, default=0.2)
    ap.add_argument("--paralogs_global", default="none",
                    help="One of: 'none', 'outliers', or a numeric threshold (e.g. 0.02).")
    ap.add_argument("--paralogs_each", choices=["yes", "no"], default="no")
    ap.add_argument("--include_hybpiper", choices=["yes","no"], default="no")
    ap.add_argument("--hybpiper_paralog_dir", "--hybpiper_dir", dest="hybpiper_paralog_dir",
                    type=Path, default=None)
    ap.add_argument("--min_len_frac_cell", type=float, default=0.0,
                    help="If > 0, require consensus_length/target_length >= this fraction per (locus,sample); "
                         "otherwise set that cell to NA before other filters.")
    ap.add_argument("--plots_pdf", action="store_true")
    ap.add_argument("--save_args_txt", action="store_true")
    args = ap.parse_args()

    subset_add = f"_{args.subset_name}" if args.subset_name else ""
    out_R = args.out_root / "00_R_objects" / args.subset_name
    out_assess = args.out_root / f"02_assessment{subset_add}"
    ensure_dir(out_R); ensure_dir(out_assess)
    if args.save_args_txt:
        (out_assess / "0_used_arguments.txt").write_text(" ".join(map(str, sys.argv)) + "\n")

    # Load from step 1a
    p_snps = out_R / "Table_SNPs.csv"
    p_len  = out_R / "Table_consensus_length.csv"
    tab_snps = read_csv_or_die(p_snps, "Table_SNPs.csv")
    tab_len  = read_csv_or_die(p_len,  "Table_consensus_length.csv")
    tab_snps.index = tab_snps.index.map(str); tab_snps.columns = tab_snps.columns.map(str)
    tab_len.index  = tab_len.index.map(str);  tab_len.columns  = tab_len.columns.map(str)

    samples = list(tab_snps.columns)

    # Targets
    target_map = parse_targets_arg(args.targets)
    max_target_length, per_genome_counts, gene_to_genome, per_genome_keyset = parse_targets_lengths(
        target_map, args.targets_format, args.target_id_mode
    )
    if max_target_length:
        print(f"[INFO] Total target loci (union): {len(max_target_length)} across {len(target_map)} FASTA(s)")
    for g, (nrec, nuniq) in per_genome_counts.items():
        print(f"[INFO] Targets parsed for {g}: records={nrec}, loci_unique={nuniq}")

    genomes = list(target_map.keys())
    genome_colors = {g:c for g,c in zip(genomes, color_cycle(len(genomes)))}

    # HybPiper cache
    hybpiper_long = []
    if args.hybpiper_paralog_dir:
        hybpiper_long = scan_hybpiper_txt(args.hybpiper_paralog_dir, samples, genomes if genomes else None)
        if hybpiper_long:
            df_hp = pd.DataFrame(hybpiper_long, columns=["sample","genome","locus"])
            df_hp.to_csv(out_R / "hybpiper_paralogs_per_sample.csv", index=False)
            print(f"[INFO] Wrote HybPiper cache: {out_R/'hybpiper_paralogs_per_sample.csv'}  (rows={len(df_hp)})")
        else:
            pd.DataFrame(columns=["sample","genome","locus"]).to_csv(out_R / "hybpiper_paralogs_per_sample.csv", index=False)
            print(f"[INFO] HybPiper TXT files not found or empty; wrote empty cache CSV.")

    # ---------- NEW: per-cell minimum length fraction filter (pre-clean) ----------
    if args.min_len_frac_cell and args.min_len_frac_cell > 0:
        common = [k for k in tab_len.index if k in max_target_length]
        tlen = pd.Series({k: max(1.0, float(max_target_length[k])) for k in common}, dtype=float)
        frac = tab_len.loc[common].div(tlen, axis=0)
        mask_low = frac < float(args.min_len_frac_cell)
        n_low = int(mask_low.values.sum())
        if n_low > 0:
            print(f"[INFO] Applying per-cell min length fraction: thr={args.min_len_frac_cell} "
                  f"| cells set to NA: {n_low}")
            tab_snps.loc[common, :] = tab_snps.loc[common, :].mask(mask_low, other=np.nan)
            tab_len.loc[common, :]  = tab_len.loc[common, :].mask(mask_low,  other=np.nan)
        else:
            print(f"[INFO] Per-cell min length fraction thr={args.min_len_frac_cell} resulted in 0 cells masked.")

    # Missing data (combined)
    loci_df = tab_snps.T
    nloci = loci_df.shape[1]; nsamples = loci_df.shape[0]

    snps_all_na_samples = set(loci_df.index[loci_df.isna().all(axis=1)].tolist())
    len_all_na_or_zero = set(
        tab_len.columns[ tab_len.isna().all(axis=0) | (tab_len.sum(axis=0, skipna=True) == 0) ].tolist()
    )
    failed_samples = sorted(snps_all_na_samples & len_all_na_or_zero)
    failed_loci = loci_df.columns[loci_df.isna().all(axis=0)]

    seq_per_locus = (~loci_df.isna()).sum(axis=0)
    seq_per_locus_prop = seq_per_locus / nsamples
    seq_per_sample = (~tab_snps.isna()).sum(axis=0)
    seq_per_sample_prop = seq_per_sample / nloci

    if max_target_length:
        loci_in_targets = [g for g in tab_len.index if g in max_target_length]
        if loci_in_targets:
            comb_target_length = sum(max_target_length[g] for g in loci_in_targets)
            comb_seq_len_samples = tab_len.loc[loci_in_targets].sum(axis=0, skipna=True)
            prop_target_length_per_sample = comb_seq_len_samples / comb_target_length

            mean_seq_length_loci = tab_len.mean(axis=1, skipna=True)
            mtl = pd.Series(np.nan, index=tab_len.index, dtype=float)
            for g in loci_in_targets:
                mtl[g] = mean_seq_length_loci[g] / float(max_target_length[g])
            prop_target_length_per_locus = mtl
        else:
            prop_target_length_per_sample = pd.Series(np.nan, index=tab_len.columns, dtype=float)
            prop_target_length_per_locus  = pd.Series(np.nan, index=tab_len.index, dtype=float)
    else:
        prop_target_length_per_sample = pd.Series(np.nan, index=tab_len.columns, dtype=float)
        prop_target_length_per_locus  = pd.Series(np.nan, index=tab_len.index, dtype=float)

    outsamples_missing_loci = seq_per_sample_prop[seq_per_sample_prop < args.thr_samples_prop_loci]
    outsamples_missing_target = prop_target_length_per_sample[prop_target_length_per_sample < args.thr_samples_prop_target]
    outsamples_missing = sorted(set(map(str, outsamples_missing_loci.index)) | set(map(str, outsamples_missing_target.index)))

    outloci_missing_samples = seq_per_locus_prop[seq_per_locus_prop < args.thr_loci_prop_samples]
    outloci_missing_target = prop_target_length_per_locus[prop_target_length_per_locus < args.thr_loci_prop_target].dropna()
    outloci_missing = sorted(set(map(str, outloci_missing_samples.index)) | set(map(str, outloci_missing_target.index)))

    # Remove low-coverage samples/loci
    tab_snps_cl1 = tab_snps.copy()
    if outsamples_missing:
        tab_snps_cl1 = tab_snps_cl1.drop(columns=outsamples_missing, errors="ignore")
    if outloci_missing:
        tab_snps_cl1 = tab_snps_cl1.drop(index=outloci_missing, errors="ignore")
    loci_cl1 = tab_snps_cl1.T

    # Map locus to genome
    locus_genome = pd.Series(index=tab_snps.index, dtype="string")
    for loc in tab_snps.index:
        locus_genome.loc[loc] = gene_to_genome.get(loc, pd.NA)

    # Plots (scatter edgecolors removed)
    def plot_missing_overview_by_genome():
        fig = plt.gcf(); fig.subplots_adjust(hspace=0.35, wspace=0.3)
        ax = plt.subplot(3, 3, 1)
        ax.boxplot(seq_per_sample_prop.dropna(), medianprops=dict(color="black", linewidth=2.5))
        ax.set_title(f"(Combined) Samples: prop. of {nloci} loci recovered")
        ax.set_xlabel(f"mean: {seq_per_sample_prop.mean():.2f} | median: {seq_per_sample_prop.median():.2f} "
                      f"| thr: {args.thr_samples_prop_loci} ({len(outsamples_missing_loci)} out)")
        ax.axhline(args.thr_samples_prop_loci, color="red", ls="--"); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 2)
        ax.boxplot(prop_target_length_per_sample.dropna(), medianprops=dict(color="black", linewidth=2.5))
        ax.set_title("(Combined) Samples: prop. of target sequence length")
        ax.set_xlabel(f"mean: {prop_target_length_per_sample.mean():.2f} | median: {prop_target_length_per_sample.median():.2f} "
                      f"| thr: {args.thr_samples_prop_target} ({len(outsamples_missing_target)} out)")
        ax.axhline(args.thr_samples_prop_target, color="red", ls="--"); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 3)
        ax.scatter(prop_target_length_per_sample, seq_per_sample_prop,
                   s=SCATTER_S_SMALL, c="black", alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
        ax.set_title("(Combined) Prop. of loci vs prop. of target length")
        ax.set_xlabel("Prop. of target length"); ax.set_ylabel("Prop. of loci")
        ax.axhline(args.thr_samples_prop_loci, color="red", ls="--")
        ax.axvline(args.thr_samples_prop_target, color="red", ls="--")
        ax.set_xlim(0,1); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 4)
        box_data = []; labels = []
        for g in genomes:
            keys = per_genome_keyset.get(g, set())
            if not keys: continue
            sub = tab_snps.loc[tab_snps.index.isin(keys)]
            if sub.shape[0] == 0: continue
            prop = (~sub.isna()).sum(axis=0) / float(sub.shape[0])
            box_data.append(prop.values); labels.append(g)
        if box_data:
            b = ax.boxplot(box_data, patch_artist=True, labels=labels,
                           medianprops=dict(color="black", linewidth=2.5))
            for patch, lab in zip(b['boxes'], labels):
                patch.set_facecolor(genome_colors[lab])
            ax.set_title("Samples: prop. of loci recovered (by genome)"); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 5)
        box_data = []; labels = []
        for g in genomes:
            keys = per_genome_keyset.get(g, set())
            keys = [k for k in keys if k in tab_len.index and k in max_target_length]
            if not keys: continue
            target_total = sum(max_target_length[k] for k in keys)
            prop = tab_len.loc[keys].sum(axis=0, skipna=True) / float(target_total if target_total>0 else np.nan)
            box_data.append(prop.values); labels.append(g)
        if box_data:
            b = ax.boxplot(box_data, patch_artist=True, labels=labels,
                           medianprops=dict(color="black", linewidth=2.5))
            for patch, lab in zip(b['boxes'], labels):
                patch.set_facecolor(genome_colors[lab])
            ax.set_title("Samples: prop. of target length (by genome)"); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 6)
        for g in genomes:
            keys = per_genome_keyset.get(g, set())
            keys = [k for k in keys if k in tab_len.index and k in tab_snps.index and k in max_target_length]
            if not keys: continue
            prop_loci = (~tab_snps.loc[keys].isna()).sum(axis=0) / float(len(keys))
            target_total = sum(max_target_length[k] for k in keys)
            prop_target = tab_len.loc[keys].sum(axis=0, skipna=True) / float(target_total if target_total>0 else np.nan)
            plt.scatter(prop_target, prop_loci, s=SCATTER_S_SMALL, label=g,
                        color=genome_colors[g], alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
        plt.title("Prop. of loci vs prop. of target length (by genome)")
        plt.xlabel("Prop. of target length"); plt.ylabel("Prop. of loci")
        plt.xlim(0,1); plt.ylim(0,1); 
        if genomes: plt.legend(fontsize=8, frameon=False)

        ax = plt.subplot(3, 3, 7)
        box_data = []; labels=[]
        for g in genomes:
            keys = [k for k in tab_snps.index if locus_genome.get(k, pd.NA) == g]
            if not keys: continue
            prop = (~tab_snps.loc[keys].T.isna()).sum(axis=0) / float(nsamples)
            box_data.append(prop.values); labels.append(g)
        if box_data:
            b = ax.boxplot(box_data, patch_artist=True, labels=labels,
                           medianprops=dict(color="black", linewidth=2.5))
            for patch, lab in zip(b['boxes'], labels):
                patch.set_facecolor(genome_colors[lab])
            ax.set_title(f"Loci: prop. of {nsamples} samples recovered (by genome)"); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 8)
        box_data = []; labels=[]
        mean_seq_length_loci = tab_len.mean(axis=1, skipna=True)
        for g in genomes:
            keys = [k for k in tab_len.index if locus_genome.get(k, pd.NA) == g and k in max_target_length]
            if not keys: continue
            vals = [mean_seq_length_loci[k] / float(max_target_length[k]) if max_target_length[k] > 0 else np.nan for k in keys]
            vals = [v for v in vals if v==v]
            if not vals: continue
            box_data.append(vals); labels.append(g)
        if box_data:
            b = ax.boxplot(box_data, patch_artist=True, labels=labels,
                           medianprops=dict(color="black", linewidth=2.5))
            for patch, lab in zip(b['boxes'], labels):
                patch.set_facecolor(genome_colors[lab])
            ax.set_title("Loci: prop. of target length (by genome)"); ax.set_ylim(0,1)

        ax = plt.subplot(3, 3, 9)
        for g in genomes:
            keys = [k for k in tab_len.index if locus_genome.get(k, pd.NA) == g and k in max_target_length]
            if not keys: continue
            prop_samples = (~tab_snps.loc[keys].T.isna()).sum(axis=0) / float(nsamples)
            prop_target = []
            for k in keys:
                v = (tab_len.loc[k].mean(skipna=True) / float(max_target_length[k])) if max_target_length[k] > 0 else np.nan
                prop_target.append(v)
            prop_target = pd.Series(prop_target, index=keys, dtype=float)
            plt.scatter(prop_target, prop_samples.loc[keys], s=SCATTER_S_SMALL,
                        color=genome_colors[g], label=g, alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
        plt.title("Prop. of samples vs prop. of target length (by genome)")
        plt.xlabel("Prop. of target length"); plt.ylabel("Prop. of samples")
        plt.xlim(0,1); plt.ylim(0,1)
        if genomes: plt.legend(fontsize=8, frameon=False)

    maybe_pdf_png(out_assess / "1_Data_recovered_overview", plot_missing_overview_by_genome, size=(15, 12), make_pdf=args.plots_pdf)

    # Per-sample/locus tables
    tab_seq_per_sample = pd.DataFrame({
        "Sample": tab_snps.columns,
        "No. loci": seq_per_sample.values,
        "Prop. of loci": seq_per_sample_prop.round(3).values,
        "Prop. of target length": prop_target_length_per_sample.round(3).values
    }).set_index("Sample")

    for g in genomes:
        keys = per_genome_keyset.get(g, set())
        keys = [k for k in keys if k in tab_snps.index and k in tab_len.index and k in max_target_length]
        if not keys:
            tab_seq_per_sample[[f"No. loci ({g})", f"Prop. of loci ({g})", f"Prop. of target length ({g})"]] = np.nan
            continue
        no_loci_g = (~tab_snps.loc[keys].isna()).sum(axis=0)
        prop_loci_g = no_loci_g / float(len(keys))
        target_total = sum(max_target_length[k] for k in keys)
        prop_tlen_g = tab_len.loc[keys].sum(axis=0, skipna=True) / float(target_total if target_total>0 else np.nan)
        tab_seq_per_sample[f"No. loci ({g})"] = no_loci_g
        tab_seq_per_sample[f"Prop. of loci ({g})"] = prop_loci_g.round(3)
        tab_seq_per_sample[f"Prop. of target length ({g})"] = prop_tlen_g.round(3)

    tab_seq_per_sample.to_csv(out_assess / "1_Data_recovered_per_sample.csv")

    per_locus = pd.DataFrame({
        "Genome": [gene_to_genome.get(l, pd.NA) for l in tab_snps.index],
        "Locus": tab_snps.index,
        "No. samples": seq_per_locus.values,
        "Prop.of samples": seq_per_locus_prop.round(3).values,
        "Prop. of target length": prop_target_length_per_locus.reindex(tab_snps.index).round(3).values
    })
    per_locus = per_locus.sort_values(["Genome", "Locus"], kind="stable").set_index("Locus")
    per_locus.to_csv(out_assess / "1_Data_recovered_per_locus.csv")

    with open(out_assess / "1_Summary_missing_data.txt", "w") as fh:
        if args.min_len_frac_cell and args.min_len_frac_cell > 0:
            fh.write(f"Per-cell min length fraction applied before filtering: {args.min_len_frac_cell}\n\n")
        fh.write("Dataset optimisation: Samples and loci removed to reduce missing data\n")
        fh.write(f"\n{len(failed_samples)} samples failed completely (SNPs all NA AND lengths all NA or zero):\n{' '.join(map(str, failed_samples))}\n")
        fh.write(f"\n{len(outsamples_missing_loci)} samples are below the threshold ({args.thr_samples_prop_loci}) for proportion of recovered loci:\n")
        for n, v in outsamples_missing_loci.items(): fh.write(f"{n}\t{v:.3f}\n")
        fh.write(f"\n{len(outsamples_missing_target)} samples are below the threshold ({args.thr_samples_prop_target}) for recovered target sequence length\n")
        for n, v in outsamples_missing_target.items(): fh.write(f"{n}\t{v:.3f}\n")
        fh.write(f"\nIn total {len(outsamples_missing)} samples were removed:\n" + "\n".join(map(str, outsamples_missing)) + "\n")
        fh.write(f"\n{len(failed_loci)} loci failed completely (SNPs all NA):\n{' '.join(map(str, failed_loci))}\n")
        fh.write(f"\n{len(outloci_missing_samples)} loci are below the threshold ({args.thr_loci_prop_samples}) for proportion of recovered samples:\n")
        for n, v in outloci_missing_samples.items(): fh.write(f"{n}\t{v:.3f}\n")
        fh.write(f"\n{len(outloci_missing_target)} loci are below the threshold ({args.thr_loci_prop_target}) for proportion of recovered target sequence length:\n")
        for n, v in outloci_missing_target.items(): fh.write(f"{n}\t{v:.3f}\n")
        fh.write(f"\nIn total {len(outloci_missing)} loci were removed:\n" + "\n".join(map(str, outloci_missing)) + "\n")

    # Step 2: paralogs (global + per-sample)
    loci_cl1 = tab_snps_cl1.T
    loci_cl1_colmeans = loci_cl1.mean(axis=0, skipna=True)
    pg = str(args.paralogs_global).strip().lower()
    if pg in {"outlier", "outliers"}:
        threshold_value = iqr_threshold(loci_cl1_colmeans.values)
        mask = loci_cl1_colmeans > threshold_value
        outloci_para_all_values = loci_cl1_colmeans[mask]
        outloci_para_all = list(outloci_para_all_values.index)
    elif pg == "none":
        threshold_value = np.inf; outloci_para_all = []; outloci_para_all_values = pd.Series(dtype=float)
    else:
        try:
            thr = float(args.paralogs_global)
        except ValueError:
            sys.exit(f"[FATAL] --paralogs_global must be 'none', 'outliers' or a number; got {args.paralogs_global}")
        threshold_value = thr
        mask = loci_cl1_colmeans > threshold_value
        outloci_para_all_values = loci_cl1_colmeans[mask]
        outloci_para_all = list(outloci_para_all_values.index)

    def plot_paralogs_all():
        fig = plt.gcf()
        ax1 = plt.subplot(2, 1, 1)
        ordered = loci_cl1_colmeans.sort_values()
        x = np.arange(len(ordered))
        colors = [genome_colors.get(gene_to_genome.get(loc, None), "gray") for loc in ordered.index]
        ax1.bar(x, ordered.values, color=colors, edgecolor="none")
        if np.isfinite(threshold_value): ax1.axhline(threshold_value, color="red", ls="--")
        removed_set = set(outloci_para_all)
        rem_x = [i for i, loc in enumerate(ordered.index) if loc in removed_set]
        rem_y = [ordered.loc[loc] for loc in ordered.index if loc in removed_set]
        if rem_x:
            ax1.scatter(rem_x, rem_y, s=16, facecolors="none", edgecolors="red", linewidths=1.4)
        ax1.set_title("Mean % SNPs across samples per locus (post missing-data filter)")
        ax1.set_xticks([])

        ax2 = plt.subplot(2, 1, 2)
        box_data = []; labels=[]
        for g in genomes:
            keys = [k for k in ordered.index if gene_to_genome.get(k, None)==g]
            if not keys: continue
            box_data.append(ordered.loc[keys].values); labels.append(g)
        if box_data:
            b = ax2.boxplot(box_data, patch_artist=True, labels=labels,
                            medianprops=dict(color="black", linewidth=2.5))
            for patch, lab in zip(b['boxes'], labels):
                patch.set_facecolor(genome_colors[lab])
        ax2.set_title("Distribution of mean % SNPs across loci (by genome)")

    maybe_pdf_png(out_assess / "2a_Paralogs_for_all_samples", plot_paralogs_all, size=(11, 7), make_pdf=args.plots_pdf)

    # Remove global paralogs
    tab_snps_cl2a = tab_snps_cl1 if len(outloci_para_all)==0 else tab_snps_cl1.drop(index=outloci_para_all, errors="ignore")

    # Per-sample (HybPhaser)
    tab_snps_cl2b = tab_snps_cl2a.copy()
    outloci_para_each = {s: [] for s in tab_snps_cl2a.columns}
    thresholds_each = {}
    tab_snps_cl2a_nozero = tab_snps_cl2a.replace(0, np.nan)
    tab_snps_cl2b_nozero = tab_snps_cl2a_nozero.copy()

    if args.paralogs_each == "yes":
        for s in tab_snps_cl2a.columns:
            series = tab_snps_cl2a_nozero[s].dropna()
            thr_s = iqr_threshold(series.values)
            thresholds_each[s] = thr_s
            removed = series[series > thr_s].index.tolist()
            outloci_para_each[s] = removed
            tab_snps_cl2b.loc[removed, s] = np.nan
            tab_snps_cl2b_nozero.loc[removed, s] = np.nan

    # HybPiper per-sample flagged set
    hybpiper_by_sample = defaultdict(set)
    if args.hybpiper_paralog_dir:
        for s, g, loc in hybpiper_long:
            hybpiper_by_sample[s].add(loc)

    def plot_paralogs_each():
        nz = tab_snps_cl2a_nozero
        order = nz.mean(axis=0, skipna=True).sort_values().index
        plt.boxplot([nz[s].dropna().values for s in order],
                    vert=False, labels=order, patch_artist=True,
                    medianprops=dict(color="black", linewidth=2.5))
        plt.xlabel("Proportion of SNPs")
        plt.title("Proportions of SNPs for all loci per sample (only loci with any SNPs)")

    maybe_pdf_png(out_assess / "2b_Paralogs_for_each_sample", plot_paralogs_each, size=(8, 14), make_pdf=args.plots_pdf)

    # Lengths aligned to cleaned SNPs
    tab_length_cl2b = tab_len.loc[tab_snps_cl2b.index, tab_snps_cl2b.columns]

    # Summary texts for paralogs
    with open(out_assess / "2_Summary_Paralogs.txt", "w") as fh:
        fh.write("Removal of putative paralog loci.\nParalogs removed for all samples:\n")
        fh.write(f"Variable 'remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs' set to: {args.paralogs_global}\n")
        if str(args.paralogs_global).lower() in {"none"}:
            fh.write("None!\n")
        else:
            if np.isfinite(threshold_value):
                fh.write(f"Resulting threshold value (mean proportion of SNPs): {threshold_value:.5f}\n")
            if len(outloci_para_all) > 0:
                fh.write(f"{len(outloci_para_all)} loci were removed:\nGenome\tLocus\tmean_prop_SNPs\n")
                for loc in order_loci_by_genome_then_name(outloci_para_all, gene_to_genome):
                    gsrc = gene_to_genome.get(loc, "")
                    fh.write(f"{gsrc}\t{loc}\t{float(outloci_para_all_values[loc]):.4f}\n")
                fh.write("\n")
        fh.write("Paralogs removed for each sample (HybPhaser vs HybPiper):\n")
        fh.write("Sample\tthreshold\t#flagged_hybphaser\t#flagged_hybpiper\t#intersection\tnames\n")
        for s in tab_snps_cl2a.columns:
            hp_set = set(outloci_para_each.get(s, []))
            pip_set = set(hybpiper_by_sample.get(s, set())) if args.include_hybpiper == "yes" else set()
            inter = hp_set & pip_set
            union_names = sorted(hp_set | pip_set)
            thr = thresholds_each.get(s, np.nan)
            fh.write(f"{s}\t{thr:.5f}\t{len(hp_set)}\t{len(pip_set)}\t{len(inter)}\t{', '.join(union_names)}\n")

    with open(out_assess / "2a_List_of_paralogs_removed_for_all_samples_with_genome.tsv", "w") as fh:
        fh.write("Genome\tLocus\tmean_prop_SNPs\n")
        for loc in order_loci_by_genome_then_name(outloci_para_all, gene_to_genome):
            gsrc = gene_to_genome.get(loc, "")
            val = float(outloci_para_all_values.get(loc, np.nan)) if len(outloci_para_all)>0 else np.nan
            fh.write(f"{gsrc}\t{loc}\t{val:.5f}\n")

    # Save cleaned core tables for 1c & downstream (sorted by Genome → Locus)
    sorted_loci = order_loci_by_genome_then_name(tab_snps_cl2b.index, gene_to_genome)
    tab_snps_cl2b_sorted = tab_snps_cl2b.loc[sorted_loci]
    tab_length_cl2b_sorted = tab_length_cl2b.loc[sorted_loci]

    tab_snps_cl2b_sorted.to_csv(out_assess / "0_Table_SNPs.csv")
    tab_length_cl2b_sorted.to_csv(out_assess / "0_Table_consensus_length.csv")
    included = [s for s in tab_snps.columns if s not in set(map(str, outsamples_missing))]
    (out_assess / "0_namelist_included_samples.txt").write_text("\n".join(included) + "\n")

    tab_snps_cl2b.to_csv(out_R / "Table_SNPs_cleaned.csv")
    tab_length_cl2b.to_csv(out_R / "Table_consensus_length_cleaned.csv")
    pd.Series(sorted(map(str, outloci_missing)), dtype="string").to_csv(out_R / "outloci_missing.csv", index=False, header=False)
    pd.Series(sorted(map(str, outsamples_missing)), dtype="string").to_csv(out_R / "outsamples_missing.csv", index=False, header=False)
    pd.Series(sorted(map(str, outloci_para_all)), dtype="string").to_csv(out_R / "outloci_para_all.csv", index=False, header=False)
    pd.DataFrame([(s, loc) for s, loci in outloci_para_each.items() for loc in loci],
                 columns=["sample", "locus"]).to_csv(out_R / "outloci_para_each.csv", index=False)

    # ---- LH & AD summary (correct n_loci_* math) ----
    if max_target_length:
        loci_kept = [g for g in tab_snps_cl2b.index if g in max_target_length]
        targets_length_cl2b = sum(max_target_length[g] for g in loci_kept) if loci_kept else np.nan
    else:
        targets_length_cl2b = np.nan

    targets_length_by_genome = {}
    for g in genomes:
        keys = [k for k in tab_snps_cl2b.index if gene_to_genome.get(k, pd.NA)==g and k in max_target_length]
        targets_length_by_genome[g] = sum(max_target_length[k] for k in keys) if keys else np.nan

    samples_final = list(tab_snps_cl2b.columns)

    base_cols = ["sample", "bp", "bpoftarget",
                 "n_loci_raw",
                 "n_loci_flagged_hybphaser_all",
                 "n_loci_flagged_hybphaser_each",
                 "n_loci_flagged_hybpiper",
                 "nloci_remaining",
                 "allele_divergence", "locus_heterozygosity",
                 "loci >0.5% SNPs", "loci >1% SNPs", "loci >2% SNPs"]
    tab_het_ad = pd.DataFrame(columns=base_cols)

    snps_pre  = tab_snps          # pre-clean (after cell-level masking if used)
    snps_post = tab_snps_cl2b     # post-clean (HybPhaser only)
    lens      = tab_length_cl2b

    for s in samples_final:
        present_pre = set(snps_pre.index[snps_pre[s].notna()])
        n_loci_raw  = len(present_pre)

        set_all  = set(outloci_para_all) & present_pre
        set_each = set(outloci_para_each.get(s, [])) & present_pre
        set_pip  = (set(hybpiper_by_sample.get(s, set())) & present_pre) if args.include_hybpiper == "yes" else set()
        flagged_union_pre = set_all | set_each | set_pip

        bp = lens[s].sum(skipna=True)
        bpoftarget = (bp / targets_length_cl2b * 100.0) if (isinstance(targets_length_cl2b,(int,float)) and targets_length_cl2b>0) else np.nan

        num = (lens[s].fillna(0) * snps_post[s].fillna(0)).sum(skipna=True)
        den = lens[s].sum(skipna=True)
        ad = (100.0 * num / den) if den > 0 else np.nan

        denom = int(snps_post[s].notna().sum())
        def lh_prop(thr):
            if denom == 0: return np.nan
            gt = (snps_post[s] > thr).fillna(False).sum()
            return round(100.0 * gt / denom, 4)

        row = {
            "sample": s,
            "bp": float(bp),
            "bpoftarget": round(bpoftarget,3) if bpoftarget==bpoftarget else np.nan,
            "n_loci_raw": n_loci_raw,
            "n_loci_flagged_hybphaser_all": len(set_all),
            "n_loci_flagged_hybphaser_each": len(set_each),
            "n_loci_flagged_hybpiper": len(set_pip),
            "nloci_remaining": n_loci_raw - len(flagged_union_pre),
            "allele_divergence": round(ad, 5) if ad==ad else np.nan,
            "locus_heterozygosity": lh_prop(0.0),
            "loci >0.5% SNPs": lh_prop(0.005),
            "loci >1% SNPs": lh_prop(0.01),
            "loci >2% SNPs": lh_prop(0.02),
        }
        tab_het_ad = pd.concat([tab_het_ad, pd.DataFrame([row])], ignore_index=True)

    # Per-genome (use pre-clean presence, per-sample)
    for g in genomes:
        keys_g = [k for k in tab_snps.index if gene_to_genome.get(k, pd.NA) == g]
        target_len_g = targets_length_by_genome.get(g, np.nan)

        col_bp = []; col_bpt = []
        col_all = []; col_each = []; col_pip = []
        col_raw = []; col_remain = []
        col_ad = []; col_lh0 = []; col_lh05 = []; col_lh1 = []; col_lh2 = []

        for s in samples_final:
            present_pre = set(snps_pre.index[snps_pre[s].notna()])
            present_pre_g = present_pre & set(keys_g)

            raw = len(present_pre_g)
            all_g  = (set(outloci_para_all) & present_pre_g)
            each_g = (set(outloci_para_each.get(s, [])) & present_pre_g)
            pip_g  = ((set(hybpiper_by_sample.get(s, set())) & present_pre_g) if args.include_hybpiper=="yes" else set())
            remain_g = raw - len(all_g | each_g | pip_g)

            sn_s = tab_snps_cl2b.loc[sorted(set(tab_snps_cl2b.index) & set(keys_g)), s] if keys_g else pd.Series(dtype=float)
            ln_s = tab_length_cl2b.loc[sorted(set(tab_length_cl2b.index) & set(keys_g)), s] if keys_g else pd.Series(dtype=float)
            bpv  = float(ln_s.sum(skipna=True)) if len(ln_s)>0 else 0.0
            col_bp.append(bpv)
            if isinstance(target_len_g,(int,float)) and target_len_g>0:
                col_bpt.append(round(bpv/target_len_g*100.0,3))
            else:
                col_bpt.append(np.nan)

            num = (ln_s.fillna(0) * sn_s.fillna(0)).sum(skipna=True)
            den = ln_s.sum(skipna=True)
            adg = (100.0 * num / den) if den > 0 else np.nan
            denom = int(sn_s.notna().sum())
            def lh_prop_g(thr):
                if denom == 0: return np.nan
                gt = (sn_s > thr).fillna(False).sum()
                return round(100.0 * gt / denom, 4)

            col_all.append(len(all_g)); col_each.append(len(each_g)); col_pip.append(len(pip_g))
            col_raw.append(raw); col_remain.append(remain_g)
            col_ad.append(round(adg,5) if adg==adg else np.nan)
            col_lh0.append(lh_prop_g(0.0))
            col_lh05.append(lh_prop_g(0.005))
            col_lh1.append(lh_prop_g(0.01))
            col_lh2.append(lh_prop_g(0.02))

        tab_het_ad[f"bp ({g})"] = col_bp
        tab_het_ad[f"bpoftarget ({g})"] = col_bpt
        tab_het_ad[f"n_loci_flagged_hybphaser_all ({g})"] = col_all
        tab_het_ad[f"n_loci_flagged_hybphaser_each ({g})"] = col_each
        tab_het_ad[f"n_loci_flagged_hybpiper ({g})"] = col_pip
        tab_het_ad[f"n_loci_raw ({g})"] = col_raw
        tab_het_ad[f"nloci_remaining ({g})"] = col_remain
        tab_het_ad[f"allele_divergence ({g})"] = col_ad
        tab_het_ad[f"locus_heterozygosity ({g})"] = col_lh0
        tab_het_ad[f"loci >0.5% SNPs ({g})"] = col_lh05
        tab_het_ad[f"loci >1% SNPs ({g})"] = col_lh1
        tab_het_ad[f"loci >2% SNPs ({g})"] = col_lh2

    tab_het_ad.to_csv(out_assess / "4_Summary_table.csv", index=False)
    tab_het_ad.to_csv(out_R / "Summary_table.csv", index=False)

    # Final plots
    def plot_lh_vs_ad_combined():
        finite = tab_het_ad[["allele_divergence", "locus_heterozygosity"]].dropna()
        plt.scatter(finite["allele_divergence"], finite["locus_heterozygosity"],
                    s=SCATTER_S_MED, c="black", alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
        plt.xlabel("Allele divergence [%]"); plt.ylabel("Locus heterozygosity [%]")
        plt.title("Locus heterozygosity vs allele divergence (combined)")

    maybe_pdf_png(out_assess / "3_LH_vs_AD", plot_lh_vs_ad_combined, size=(8, 8), make_pdf=args.plots_pdf)

    def plot_var_lh_grid_combined():
        fig = plt.gcf(); fig.subplots_adjust(hspace=0.35, wspace=0.3)
        panels = [("locus_heterozygosity", "LH (0% SNPs)"),
                  ("loci >0.5% SNPs", "LH (>0.5% SNPs)"),
                  ("loci >1% SNPs", "LH (>1% SNPs)"),
                  ("loci >2% SNPs", "LH (>2% SNPs)")]
        for i, (col, title) in enumerate(panels, start=1):
            ax = plt.subplot(2, 2, i)
            finite = tab_het_ad[["allele_divergence", col]].dropna()
            ax.scatter(finite["allele_divergence"], finite[col],
                       s=SCATTER_S_SMALL, c="black", alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
            ax.set_xlabel("Allele divergence [%]"); ax.set_ylabel(f"{title} [%]")
            ax.set_title(f"{title} vs allele divergence (combined)")

    maybe_pdf_png(out_assess / "3_varLH_vs_AD", plot_var_lh_grid_combined, size=(8, 8), make_pdf=args.plots_pdf)

    def plot_lh_vs_ad_by_genome():
        for g in genomes:
            x = tab_het_ad[f"allele_divergence ({g})"]
            y = tab_het_ad[f"locus_heterozygosity ({g})"]
            plt.scatter(x, y, s=SCATTER_S_MED, color=genome_colors[g], label=g, alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
        plt.xlabel("Allele divergence [%]"); plt.ylabel("Locus heterozygosity [%]")
        plt.title("Locus heterozygosity vs allele divergence (by genome)")
        if genomes: plt.legend(frameon=False, fontsize=8)

    maybe_pdf_png(out_assess / "3_LH_vs_AD_by_genome", plot_lh_vs_ad_by_genome, size=(8, 8), make_pdf=args.plots_pdf)

    def plot_var_lh_grid_by_genome():
        fig = plt.gcf(); fig.subplots_adjust(hspace=0.35, wspace=0.3)
        panels = [("locus_heterozygosity", "LH (0% SNPs)"),
                  ("loci >0.5% SNPs", "LH (>0.5% SNPs)"),
                  ("loci >1% SNPs", "LH (>1% SNPs)"),
                  ("loci >2% SNPs", "LH (>2% SNPs)")]
        for i, (base, title) in enumerate(panels, start=1):
            ax = plt.subplot(2, 2, i)
            for g in genomes:
                x = tab_het_ad[f"allele_divergence ({g})"]
                y = tab_het_ad[f"{base} ({g})"]
                ax.scatter(x, y, s=SCATTER_S_SMALL, color=genome_colors[g], label=g if i==1 else None, alpha=SCATTER_ALPHA, edgecolors='none', linewidths=0)
            ax.set_xlabel("Allele divergence [%]"); ax.set_ylabel(f"{title} [%]")
            ax.set_title(f"{title} vs allele divergence (by genome)")
            if i==1 and genomes: ax.legend(frameon=False, fontsize=8)

    maybe_pdf_png(out_assess / "3_varLH_vs_AD_by_genome", plot_var_lh_grid_by_genome, size=(8, 8), make_pdf=args.plots_pdf)

    print(f"[INFO] Done. Wrote outputs to:\n  {out_assess}\n  {out_R}")

if __name__ == "__main__":
    main()
