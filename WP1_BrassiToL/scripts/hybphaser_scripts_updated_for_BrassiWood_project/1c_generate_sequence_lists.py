#!/usr/bin/env python3
# 1c_generate_sequence_lists.py
#
# Port of HybPhaser 1c (R) to Python.
# Builds cleaned sequence lists (per locus and per sample) for consensus and contigs,
# applying removal decisions produced by 1b (missing data + paralogs).
#
# Inputs expected (from step 1b): <out_root>/00_R_objects/<subset>/
#   - Table_SNPs_cleaned.csv
#   - outloci_missing.csv
#   - outsamples_missing.csv
#   - outloci_para_all.csv
#   - outloci_para_each.csv          (columns: sample,locus)
#   - hybpiper_paralogs_per_sample.csv   (optional; columns: sample,genome,locus)
#
# FASTA inputs (flat layout):
#   consensus_dir: files like <sample>_<genome>_consensus.fasta  (multi-locus)
#   contig_dir   : files like <sample>_<genome>.fasta            (HybPiper; multi-locus)
#
# Targets (same as 1a/1b):
#   --targets GENOME=path ...  + --targets_format + --target_id_mode
#   We only keep loci present in the target files ("target loci").
#
# Outputs:
#   <out_root>/03_sequence_lists_<subset>/
#     ├─ loci_consensus/   <locus>_consensus.fasta
#     ├─ loci_contigs/     <locus>_contig.fasta
#     ├─ samples_consensus/<sample>_consensus.fasta
#     ├─ samples_contigs/  <sample>_contig.fasta
#     ├─ overview_summary.txt
#     ├─ per_sample_paralog_counts.tsv
#     ├─ loci_removed_global.txt
#     ├─ loci_written_counts.tsv
#     └─ (from earlier) 00_R_objects/<subset>/missing_per_sample_inputs_1c.txt
#
# Notes:
# - Sequence IDs in per-locus FASTAs are sample names. In per-sample FASTAs they are "<sample>-<locus>".
# - "--intronerated yes" only changes filenames to match HybPhaser naming (content unchanged).

import argparse
from pathlib import Path
import sys
import re
from collections import defaultdict, Counter
from typing import Optional, List, Dict, Set, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

# ---------- helpers ----------

FA_EXTS = (".fa", ".fas", ".fasta", ".fna")

FILE_RE_CONS = re.compile(
    r"^(?P<sample>.+)_(?P<genome>[^_]+)_consensus\.(?:fa|fas|fasta|fna)$",
    re.IGNORECASE
)
FILE_RE_CONTIG = re.compile(
    r"^(?P<sample>.+)_(?P<genome>[^_]+)\.(?:fa|fas|fasta|fna)$",
    re.IGNORECASE
)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def seq_writer_unwrapped(handle, records: List[SeqRecord]):
    w = FastaWriter(handle, wrap=None)
    w.write_file(records)

def infer_locus_from_id(rec_id: str, fallback: str) -> str:
    rid = re.sub(r"\s.*$", "", rec_id)
    if "-" in rid:
        return rid.rsplit("-", 1)[-1] or fallback
    return rid or fallback

# ---------- IO ----------

def read_index_col_csv(p: Path) -> pd.DataFrame:
    if not p.exists():
        sys.exit(f"[FATAL] Missing required CSV: {p}")
    return pd.read_csv(p, index_col=0, low_memory=False)

def read_single_col_list(p: Path) -> List[str]:
    """Read a 1-column CSV or plain-text list. Returns [] if file is missing or empty."""
    if not p.exists():
        return []
    try:
        if p.stat().st_size == 0:
            return []
    except OSError:
        return []
    # First try: read as plain text, one item per non-empty line
    try:
        lines = [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]
        if lines:
            return lines
    except Exception:
        pass
    # Fallback: CSV
    try:
        df = pd.read_csv(p, header=None)
        if df.empty:
            return []
        return [str(x) for x in df.iloc[:, 0].dropna().astype(str).tolist()]
    except pd.errors.EmptyDataError:
        return []
    except Exception:
        # Be forgiving: treat unreadable as empty list
        return []


def read_long_two_col(p: Path, col1="sample", col2="locus") -> Dict[str, Set[str]]:
    """Read a 2-column long table (CSV). Returns {} if missing/empty or columns absent."""
    if not p.exists():
        return {}
    try:
        if p.stat().st_size == 0:
            return {}
    except OSError:
        return {}
    try:
        df = pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return {}
    except Exception:
        return {}
    if df.empty or col1 not in df.columns or col2 not in df.columns:
        return {}
    out = defaultdict(set)
    for _, r in df.iterrows():
        s = str(r[col1]) if pd.notna(r[col1]) else ""
        l = str(r[col2]) if pd.notna(r[col2]) else ""
        if s and l:
            out[s].add(l)
    return dict(out)


# ---------- targets ----------

def gene_key_from_header(header: str, mode: str) -> str:
    if mode in ("dash", "full"):
        rid = re.sub(r"\s.*$", "", header)
    elif mode == "space":
        rid = header.split()[0]
    else:
        rid = re.sub(r"\s.*$", "", header)
    if mode == "dash" or mode == "space":
        return rid.rsplit("-", 1)[-1] if "-" in rid else rid
    return rid

def parse_targets_arg(items: List[str]) -> Dict[str, Path]:
    d = {}
    for it in items:
        if "=" not in it:
            sys.exit(f"[FATAL] --targets entry must be GENOME=/path/to.fasta, got '{it}'")
        k, v = it.split("=", 1)
        k = k.strip()
        v = v.strip()
        d[k] = Path(v)
    return d

def parse_targets_lengths(target_map: Dict[str, Path], fmt: str, id_mode: str):
    is_aa = fmt.upper() == "AA"
    max_per_gene: Dict[str, int] = {}
    counts: Dict[str, Tuple[int, int]] = {}
    for genome, fasta_path in target_map.items():
        nrecs = 0
        seen = {}
        if not fasta_path.exists():
            counts[genome] = (0, 0)
            continue
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            nrecs += 1
            key = gene_key_from_header(rec.id, id_mode)
            L = len(rec.seq) * (3 if is_aa else 1)
            if key not in seen or L > seen[key]:
                seen[key] = L
            if key not in max_per_gene or L > max_per_gene[key]:
                max_per_gene[key] = L
        counts[genome] = (nrecs, len(seen))
    return max_per_gene, counts

# ---------- FASTA scanning ----------

def scan_flat_fastas(root: Optional[Path], samples_set: Set[str], kind: str) -> Dict[str, List[Path]]:
    """
    Returns {sample: [Path,...]} for given kind ('consensus' with _consensus suffix,
    or 'contig' without).
    """
    bucket: Dict[str, List[Path]] = defaultdict(list)
    if not root or not root.exists():
        return {}

    for fp in root.rglob("*"):
        if not fp.is_file() or fp.suffix.lower() not in FA_EXTS:
            continue

        if kind == "consensus":
            m = FILE_RE_CONS.match(fp.name)
        else:
            # avoid misclassifying "*_consensus.fasta" as contigs
            if fp.stem.lower().endswith("_consensus"):
                continue
            m = FILE_RE_CONTIG.match(fp.name)

        if not m:
            continue

        smp = m.group("sample")
        if smp in samples_set:
            bucket[smp].append(fp)

    return dict(bucket)

def parse_per_sample_fastas(file_list: List[Path], sample: str, target_loci: Set[str]) -> Dict[str, SeqRecord]:
    loci_map: Dict[str, SeqRecord] = {}
    for fp in file_list:
        for rec in SeqIO.parse(str(fp), "fasta"):
            locus = infer_locus_from_id(rec.id, fallback=fp.stem)
            if locus not in target_loci:
                continue
            loci_map[locus] = SeqRecord(Seq(str(rec.seq)), id=sample, description="")
    return loci_map

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="HybPhaser 1c in Python (sequence list generation)")
    ap.add_argument("--out_root", required=True, type=Path)
    ap.add_argument("--subset_name", default="")
    ap.add_argument("--samples_file", required=True, type=Path)
    ap.add_argument("--intronerated", choices=["yes","no"], default="no")
    ap.add_argument("--layout", choices=["flat"], default="flat")
    ap.add_argument("--consensus_dir", required=True, type=Path)
    ap.add_argument("--contig_dir", type=Path, default=None)
    ap.add_argument("--remove_hybpiper", choices=["yes","no"], default="yes")

    # targets (filter to target loci only)
    ap.add_argument("--targets", action="append", default=[], help="Repeat: GENOME=/path/to/targets.fasta")
    ap.add_argument("--targets_format", choices=["DNA","AA"], default="DNA")
    ap.add_argument("--target_id_mode", choices=["dash","space","full"], default="dash")

    ap.add_argument("--log_every", type=int, default=200)

    args = ap.parse_args()
    subset_add = f"_{args.subset_name}" if args.subset_name else ""
    out_R = args.out_root / "00_R_objects" / args.subset_name
    out_seq = args.out_root / f"03_sequence_lists{subset_add}"

    # output dirs
    loci_consensus_dir = out_seq / "loci_consensus"
    loci_contigs_dir   = out_seq / "loci_contigs"
    samples_cons_dir   = out_seq / "samples_consensus"
    samples_cont_dir   = out_seq / "samples_contigs"
    for d in (loci_consensus_dir, loci_contigs_dir, samples_cons_dir, samples_cont_dir):
        ensure_dir(d)

    intron_name = "intronerated" if args.intronerated == "yes" else ""
    intron_us   = "_" if intron_name else ""

    # load samples
    samples = [ln.strip() for ln in args.samples_file.read_text().splitlines() if ln.strip()]
    samples_set = set(samples)

    # --- read decisions from 1b ---
    tab_snps_cl2b = read_index_col_csv(out_R / "Table_SNPs_cleaned.csv")   # loci x samples
    failed_loci = tab_snps_cl2b.index[tab_snps_cl2b.isna().all(axis=1)].map(str).tolist()

    outloci_missing    = read_single_col_list(out_R / "outloci_missing.csv")
    outsamples_missing = read_single_col_list(out_R / "outsamples_missing.csv")
    outloci_para_all   = read_single_col_list(out_R / "outloci_para_all.csv")
    para_each_hphase   = read_long_two_col(out_R / "outloci_para_each.csv")  # {sample: set(loci)}

    para_each_hpiper: Dict[str, Set[str]] = {}
    if args.remove_hybpiper == "yes":
        hp = out_R / "hybpiper_paralogs_per_sample.csv"
        if hp.exists():
            para_each_hpiper = read_long_two_col(hp)

    # per-sample removal union
    per_sample_rm: Dict[str, Set[str]] = defaultdict(set)
    for s in tab_snps_cl2b.columns:
        if s in para_each_hphase:
            per_sample_rm[s].update(para_each_hphase[s])
        if s in para_each_hpiper:
            per_sample_rm[s].update(para_each_hpiper[s])

    samples_included = [s for s in samples if s not in set(outsamples_missing)]
    loci_remove_all = set(failed_loci) | set(outloci_missing) | set(outloci_para_all)

    # --- targets: keep only target loci ---
    target_map = parse_targets_arg(args.targets)
    max_target_len, counts = parse_targets_lengths(target_map, args.targets_format, args.target_id_mode)
    target_loci = set(max_target_len.keys())  # previously named "allowed_loci"

    # --- scan input FASTAs ---
    cons_index  = scan_flat_fastas(args.consensus_dir, samples_set, "consensus")
    contig_index= scan_flat_fastas(args.contig_dir,    samples_set, "contig") if args.contig_dir else {}

    missing_inputs = []

    # parse consensus sequences
    seqs_consensus: Dict[str, Dict[str, SeqRecord]] = defaultdict(dict)  # locus -> {sample: rec}
    for s in samples_included:
        files_s = cons_index.get(s, [])
        if not files_s:
            missing_inputs.append(f"consensus:{s}")
            continue
        loci_map = parse_per_sample_fastas(files_s, s, target_loci)
        for loc, rec in loci_map.items():
            seqs_consensus[loc][s] = rec

    # parse contigs (optional)
    seqs_contigs: Dict[str, Dict[str, SeqRecord]] = defaultdict(dict)
    if args.contig_dir:
        for s in samples_included:
            files_s = contig_index.get(s, [])
            if not files_s:
                missing_inputs.append(f"contig:{s}")
                continue
            loci_map = parse_per_sample_fastas(files_s, s, target_loci)
            for loc, rec in loci_map.items():
                seqs_contigs[loc][s] = rec

    if missing_inputs:
        ensure_dir(out_R)
        (out_R / "missing_per_sample_inputs_1c.txt").write_text("\n".join(sorted(set(missing_inputs))) + "\n")

    # --- write outputs & collect stats ---
    loci_written_counts: Dict[str, int] = {}
    loci_written_counts_contig: Dict[str, int] = {}
    samples_with_consensus_written = 0
    samples_with_contigs_written   = 0

    # per-locus consensus
    n_locus_cons = 0
    for locus, smap in seqs_consensus.items():
        if locus in loci_remove_all:
            continue
        recs: List[SeqRecord] = []
        for s in samples_included:
            if s in smap and locus not in per_sample_rm[s]:
                recs.append(SeqRecord(smap[s].seq, id=s, description=""))
        if recs:
            outp = (loci_consensus_dir / f"{locus}_{intron_name}{intron_us}consensus.fasta") if intron_name else (loci_consensus_dir / f"{locus}_consensus.fasta")
            with open(outp, "w") as h:
                seq_writer_unwrapped(h, recs)
            n_locus_cons += 1
            loci_written_counts[locus] = len(recs)

    # per-locus contigs
    n_locus_contig = 0
    if seqs_contigs:
        for locus, smap in seqs_contigs.items():
            if locus in loci_remove_all:
                continue
            recs: List[SeqRecord] = []
            for s in samples_included:
                if s in smap and locus not in per_sample_rm[s]:
                    recs.append(SeqRecord(smap[s].seq, id=s, description=""))
            if recs:
                outp = (loci_contigs_dir / f"{locus}_{intron_name}{intron_us}contig.fasta") if intron_name else (loci_contigs_dir / f"{locus}_contig.fasta")
                with open(outp, "w") as h:
                    seq_writer_unwrapped(h, recs)
                n_locus_contig += 1
                loci_written_counts_contig[locus] = len(recs)

    # per-sample consensus
    n_sample_cons = 0
    for i, s in enumerate(samples_included, start=1):
        recs: List[SeqRecord] = []
        rm = per_sample_rm.get(s, set())
        for locus, smap in seqs_consensus.items():
            if locus in loci_remove_all or locus in rm:
                continue
            if s in smap:
                recs.append(SeqRecord(smap[s].seq, id=f"{s}-{locus}", description=""))
        if recs:
            outp = (samples_cons_dir / f"{s}{intron_us}{intron_name}_consensus.fasta") if intron_name else (samples_cons_dir / f"{s}_consensus.fasta")
            with open(outp, "w") as h:
                seq_writer_unwrapped(h, recs)
            n_sample_cons += 1
            samples_with_consensus_written += 1
        if (i % max(1, args.log_every) == 0) or (i == len(samples_included)):
            print(f"[INFO]   Per-sample consensus: processed {i}/{len(samples_included)}")

    # per-sample contigs
    n_sample_contig = 0
    if seqs_contigs:
        for i, s in enumerate(samples_included, start=1):
            recs: List[SeqRecord] = []
            rm = per_sample_rm.get(s, set())
            for locus, smap in seqs_contigs.items():
                if locus in loci_remove_all or locus in rm:
                    continue
                if s in smap:
                    recs.append(SeqRecord(smap[s].seq, id=f"{s}-{locus}", description=""))
            if recs:
                outp = (samples_cont_dir / f"{s}{intron_us}{intron_name}_contig.fasta") if intron_name else (samples_cont_dir / f"{s}_contig.fasta")
                with open(outp, "w") as h:
                    seq_writer_unwrapped(h, recs)
                n_sample_contig += 1
                samples_with_contigs_written += 1
            if (i % max(1, args.log_every) == 0) or (i == len(samples_included)):
                print(f"[INFO]   Per-sample contigs : processed {i}/{len(samples_included)}")

    # --- per-sample paralog counts (for summary & TSV) ---
    rows_ps = []
    total_union_set = set()
    total_hph_set = set()
    total_hpp_set = set()
    for s in samples_included:
        c_hph = len(para_each_hphase.get(s, set()))
        c_hpp = len(para_each_hpiper.get(s, set()))
        union_set = set()
        union_set.update(para_each_hphase.get(s, set()))
        union_set.update(para_each_hpiper.get(s, set()))
        c_union = len(union_set)
        rows_ps.append((s, c_hph, c_hpp, c_union))
        total_union_set.update(union_set)
        total_hph_set.update(para_each_hphase.get(s, set()))
        total_hpp_set.update(para_each_hpiper.get(s, set()))

    df_ps = pd.DataFrame(rows_ps, columns=["sample", "paralogs_each_hybphaser", "paralogs_each_hybpiper", "paralogs_each_union"])
    if not df_ps.empty:
        df_ps.sort_values(by=["paralogs_each_union","paralogs_each_hybphaser","paralogs_each_hybpiper","sample"], ascending=[False,False,False,True], inplace=True)
        df_ps.to_csv(out_seq / "per_sample_paralog_counts.tsv", sep="\t", index=False)

    # write loci removed list
    if loci_remove_all:
        (out_seq / "loci_removed_global.txt").write_text("\n".join(sorted(loci_remove_all)) + "\n")

    # write loci_written_counts
    if loci_written_counts or loci_written_counts_contig:
        df_lwc = pd.DataFrame({
            "locus": sorted(set(list(loci_written_counts.keys()) + list(loci_written_counts_contig.keys())))
        })
        df_lwc["n_samples_consensus"] = df_lwc["locus"].map(lambda x: loci_written_counts.get(x, 0))
        df_lwc["n_samples_contig"]    = df_lwc["locus"].map(lambda x: loci_written_counts_contig.get(x, 0))
        df_lwc.to_csv(out_seq / "loci_written_counts.tsv", sep="\t", index=False)

    # --- overview summary ---
    n_total_samples = len(samples)
    n_included_samples = len(samples_included)
    n_outsamples_missing = len(outsamples_missing)

    n_failed_loci = len(set(failed_loci))
    n_outloci_missing = len(set(outloci_missing))
    n_outloci_para_all = len(set(outloci_para_all))
    n_loci_remove_all = len(loci_remove_all)

    # missing inputs tallies
    miss_consensus = sum(1 for x in set(missing_inputs) if x.startswith("consensus:"))
    miss_contig    = sum(1 for x in set(missing_inputs) if x.startswith("contig:"))

    # top 10 samples by union paralogs
    top_lines = []
    if not df_ps.empty:
        top = df_ps.head(10).itertuples(index=False)
        for row in top:
            top_lines.append(f"  {row.sample}\tunion={row.paralogs_each_union}\t(hybphaser={row.paralogs_each_hybphaser}, hybpiper={row.paralogs_each_hybpiper})")

    summary_lines = []
    summary_lines.append(f"Samples included: {n_included_samples} / {n_total_samples}")
    summary_lines.append(f"Samples removed globally (missing-data thresholds from 1b): {n_outsamples_missing}")
    summary_lines.append("")
    summary_lines.append("Global locus removals (applied to all samples):")
    summary_lines.append(f"  failed_loci (all-NA after cleaning): {n_failed_loci}")
    summary_lines.append(f"  outloci_missing (missing-data thresholds): {n_outloci_missing}")
    summary_lines.append(f"  outloci_para_all (global paralogs): {n_outloci_para_all}")
    summary_lines.append(f"  union removed (failed + missing + paralogs_all): {n_loci_remove_all}")
    summary_lines.append("")
    summary_lines.append("Per-sample paralog removals (counts per sample):")
    summary_lines.append(f"  unique loci flagged by HybPhaser across samples: {len(total_hph_set)}")
    summary_lines.append(f"  unique loci flagged by HybPiper  across samples: {len(total_hpp_set)}")
    summary_lines.append(f"  unique loci flagged by union (HybPhaser ∪ HybPiper): {len(total_union_set)}")
    if top_lines:
        summary_lines.append("  Top 10 samples by union paralogs:")
        summary_lines.extend(top_lines)
    else:
        summary_lines.append("  (no per-sample paralogs found)")
    summary_lines.append("")
    summary_lines.append("Outputs written:")
    summary_lines.append(f"  per-locus consensus FASTAs: {n_locus_cons}")
    summary_lines.append(f"  per-locus contig FASTAs   : {n_locus_contig}")
    summary_lines.append(f"  per-sample consensus FASTAs: {n_sample_cons}")
    summary_lines.append(f"  per-sample contig FASTAs   : {n_sample_contig}")
    summary_lines.append("")
    summary_lines.append("Missing input files (per-sample):")
    summary_lines.append(f"  consensus missing: {miss_consensus}")
    summary_lines.append(f"  contig missing   : {miss_contig}")
    if missing_inputs:
        summary_lines.append("  (full list in 00_R_objects/<subset>/missing_per_sample_inputs_1c.txt)")

    (out_seq / "overview_summary.txt").write_text("\n".join(summary_lines) + "\n")

    print(f"[INFO] Done. Outputs written to {out_seq}")

if __name__ == "__main__":
    main()
