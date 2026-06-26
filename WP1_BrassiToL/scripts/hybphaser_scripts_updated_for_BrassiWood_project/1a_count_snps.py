#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HybPhaser 1a (SNP counting) in Python — flat layout + multi-genome.

Behavior (robust filenames):
- For each sample in --samples_file and each genome in your targets mapping,
  we look for (in this priority order):
    1) <sample>_<genome>_consensus.(fa|fas|fasta|fna)
    2) <sample>_<genome>.(fa|fas|fasta|fna)
- Search is recursive within --consensus_dir.
- Sample IDs may contain underscores; this logic does not use regex splitting.

Examples accepted:
  SRR8528391_B764_consensus.fasta
  S_670_B764_consensus.fasta
  23Q02_A353.fasta

Outputs:
  out_root/00_R_objects/<subset_name or combined>/
      Table_SNPs.csv
      Table_consensus_length.csv
      targets_file_summary.csv
      consensus_files_per_sample.csv
      locus_source_map.csv
      missing_consensus_files.txt          (expected but not found)
      duplicate_stems_in_consensus_dir.txt (diagnostic; stems seen >1 time)
"""

from __future__ import annotations
import argparse, sys
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from Bio import SeqIO

# --- constants ---
FA_EXTS = {".fa", ".fas", ".fasta", ".fna"}
AMBIG_1 = set("YKRSMWykrsmaw")
AMBIG_2 = set("DHVBdhvb")
STRIP   = set("Nn?-")


def parse_args():
    p = argparse.ArgumentParser(
        description="Count SNP proportions (ambiguity codes) per locus per sample, across one or many genomes."
    )
    p.add_argument("--consensus_dir", required=True, type=Path,
                   help="Folder with <sample>_<genome>[_consensus].fasta files (recursive).")
    p.add_argument("--samples_file", required=True, type=Path,
                   help="Text file with one sample per line (exactly as used in filenames).")
    # Target FASTAs: either explicit mapping or legacy paired lists
    p.add_argument("--targets", action="append", default=[],
                   help="Repeatable: GENOME=/path/to/targets.fasta (recommended).")
    p.add_argument("--genomes", nargs="*", default=[],
                   help="Legacy mode: list of genomes (must match order/length of --targets_fasta).")
    p.add_argument("--targets_fasta", nargs="*", default=[],
                   help="Legacy mode: list of target FASTAs aligned to --genomes by position.")
    p.add_argument("--targets_format", choices=["DNA","AA"], default="DNA",
                   help="Targets FASTA format; AA lengths ×3 (only affects length summaries downstream).")
    p.add_argument("--target_id_mode", choices=["dash","full"], default="dash",
                   help="dash: use last token after '-' ; full: use entire header token up to whitespace.")
    p.add_argument("--allow_extra_loci", choices=["yes","no"], default="no",
                   help="Include loci present in consensus but NOT in targets [default: no].")
    p.add_argument("--out_root", required=True, type=Path, help="Output root directory.")
    p.add_argument("--subset_name", default="", help="Optional label for 00_R_objects/<subset_name>. Otherwise 'combined'.")
    p.add_argument("--log_every", type=int, default=200, help="Print progress every N files [default 200].")
    return p.parse_args()


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def infer_locus(header: str, mode: str) -> str:
    """Extract a locus id from a FASTA header depending on target_id_mode."""
    # take the first whitespace-separated token
    rid = header.split()[0]
    if mode == "full":
        return rid
    # "dash": last token after '-'
    return rid.split("-")[-1] if "-" in rid else rid


def clean_len_and_prop(seq: str) -> tuple[int, float | float]:
    # length excludes N/n/?/-
    raw = [c for c in str(seq) if c not in STRIP]
    L = len(raw)
    if L == 0:
        return 0, np.nan
    c1 = sum(1 for c in raw if c in AMBIG_1)
    c2 = sum(1 for c in raw if c in AMBIG_2)
    ambigs = c1 + 2 * c2
    return L, (ambigs / L)


def parse_targets_mapping(args) -> tuple[dict[str, Path], dict[str, set[str]], dict[str, str]]:
    """
    Returns:
      tmap : genome -> Path(fasta)
      loci_per_genome : genome -> set(loci)
      locus_source : locus -> genome (first genome to define that locus; avoids duplicates)
    """
    tmap: dict[str, Path] = {}
    loci_per_genome: dict[str, set[str]] = {}
    locus_source: dict[str, str] = {}

    # explicit mapping takes precedence
    if args.targets:
        for item in args.targets:
            if "=" not in item:
                sys.exit(f"[FATAL] --targets must be GENOME=/path/to/file.fasta ; got: {item}")
            genome, path = item.split("=", 1)
            genome = genome.strip()
            if not genome:
                sys.exit("[FATAL] empty genome in --targets")
            p = Path(path).expanduser().resolve()
            if not p.exists():
                sys.exit(f"[FATAL] targets FASTA not found: {p}")
            tmap[genome] = p
    else:
        # legacy paired lists
        if not args.genomes or not args.targets_fasta or (len(args.genomes) != len(args.targets_fasta)):
            sys.exit("[FATAL] Use either --targets GENOME=FASTA... OR paired --genomes ... --targets_fasta ... (same length).")
        for g, f in zip(args.genomes, args.targets_fasta):
            p = Path(f).expanduser().resolve()
            if not p.exists():
                sys.exit(f"[FATAL] targets FASTA not found: {p}")
            tmap[g] = p

    # parse targets
    print(f"[INFO] Parsing target FASTAs ({len(tmap)}) ...")
    for i, (genome, fasta) in enumerate(tmap.items(), start=1):
        print(f"[INFO]   [{i}/{len(tmap)}] {genome}: {fasta}")
        loci = set()
        nrec = 0
        try:
            for rec in SeqIO.parse(str(fasta), "fasta"):
                nrec += 1
                locus = infer_locus(rec.id, args.target_id_mode)
                loci.add(locus)
                locus_source.setdefault(locus, genome)
        except Exception as e:
            sys.exit(f"[FATAL] Failed reading targets FASTA {fasta}: {e}")
        print(f"[INFO]        records={nrec}, loci_parsed={len(loci)}")
        if len(loci) == 0:
            sys.exit(f"[FATAL] Parsed 0 loci from targets FASTA: {fasta}\n       Try --target_id_mode full if headers differ.")
        loci_per_genome[genome] = loci

    union = set().union(*loci_per_genome.values()) if loci_per_genome else set()
    print(f"[INFO] Total target loci (union): {len(union)}")
    return tmap, loci_per_genome, locus_source


def build_fasta_index(consensus_dir: Path):
    """
    Walk once; build:
      - index: stem_lower -> [Path, ...] (stem = filename without extension)
      - stats: counts and duplicate stems for diagnostics
    """
    index: dict[str, list[Path]] = defaultdict(list)
    n_paths = 0
    for fp in consensus_dir.rglob("*"):
        if fp.is_file():
            ext = fp.suffix.lower()
            if ext in FA_EXTS:
                n_paths += 1
                index[fp.stem.lower()].append(fp)

    dup_stems = {stem: paths for stem, paths in index.items() if len(paths) > 1}
    print(f"[INFO] Indexed {n_paths} FASTA files under {consensus_dir} (unique stems: {len(index)})")
    if dup_stems:
        print(f"[WARN] Found {len(dup_stems)} duplicate stems (same name, multiple locations). See diagnostic file.")
    return index, dup_stems


def resolve_consensus_path(index: dict[str, list[Path]], sample: str, genome: str) -> Path | None:
    """
    Choose file by priority:
      1) <sample>_<genome>_consensus
      2) <sample>_<genome>
    """
    cand1 = f"{sample}_{genome}_consensus".lower()
    cand2 = f"{sample}_{genome}".lower()
    for stem in (cand1, cand2):
        paths = index.get(stem, [])
        if paths:
            # If multiple, pick the first deterministically sorted by path
            return sorted(paths)[0]
    return None


def main():
    args = parse_args()

    out_sub = args.subset_name if args.subset_name else "combined"
    out_R = args.out_root / "00_R_objects" / out_sub
    ensure_dir(out_R)

    # Read samples
    samples = [ln.strip() for ln in args.samples_file.read_text().splitlines() if ln.strip()]
    if not samples:
        sys.exit("[FATAL] samples_file is empty")
    print(f"[INFO] Loaded {len(samples)} sample IDs from {args.samples_file}")

    # Targets
    tmap, loci_per_genome, locus_source = parse_targets_mapping(args)
    target_union = set().union(*loci_per_genome.values()) if loci_per_genome else set()

    # Save target summary
    pd.DataFrame(
        [{"genome": g, "targets_fasta": str(p), "n_loci": len(loci_per_genome[g])} for g, p in tmap.items()]
    ).to_csv(out_R / "targets_file_summary.csv", index=False)
    pd.DataFrame(
        [{"locus": l, "genome": locus_source.get(l, "")} for l in sorted(target_union)]
    ).to_csv(out_R / "locus_source_map.csv", index=False)

    # Build FASTA index once
    index, dup_stems = build_fasta_index(args.consensus_dir)
    # Write duplicate stems diagnostic (if any)
    if dup_stems:
        with open(out_R / "duplicate_stems_in_consensus_dir.txt", "w") as fh:
            for stem, paths in sorted(dup_stems.items()):
                for p in sorted(paths):
                    fh.write(f"{stem}\t{p}\n")

    # Resolve expected files deterministically
    genomes_filter = list(tmap.keys())
    files = []
    missing = []
    for s in samples:
        for g in genomes_filter:
            path = resolve_consensus_path(index, s, g)
            if path is None:
                # Record both expected stems (so it’s crystal clear)
                missing.append(f"{s}_{g}_consensus.(fa|fas|fasta|fna)")
                missing.append(f"{s}_{g}.(fa|fas|fasta|fna)")
            else:
                files.append((path, s, g))

    # Save missing list (still proceed with found files)
    if missing:
        miss_file = out_R / "missing_consensus_files.txt"
        with open(miss_file, "w") as fh:
            for line in sorted(set(missing)):
                fh.write(line + "\n")
        print(f"[WARN] {len(set(missing))} expected consensus filenames were not found. See: {miss_file}")

    if not files:
        sys.exit("[FATAL] No consensus files resolved. Check naming and --consensus_dir.")

    print(f"[INFO] Using {len(files)} consensus files (resolved from sample×genome pairs).")

    # Build sparse maps: locus -> sample -> metric
    snp_map: dict[str, dict[str, float]] = defaultdict(dict)
    len_map: dict[str, dict[str, float]] = defaultdict(dict)

    # Parse per file
    for i, (fp, smp, gen) in enumerate(files, start=1):
        if i % max(1, args.log_every) == 0 or i == len(files):
            print(f"[INFO]   parsing {i}/{len(files)}: {fp.name}")
        try:
            for rec in SeqIO.parse(str(fp), "fasta"):
                locus = infer_locus(rec.id, args.target_id_mode)
                # skip non-target loci unless allowed
                if locus not in target_union and args.allow_extra_loci == "no":
                    continue
                L, prop = clean_len_and_prop(rec.seq)
                len_map[locus][smp] = L
                snp_map[locus][smp] = prop
        except Exception as e:
            print(f"[WARN] Failed to parse {fp}: {e}", file=sys.stderr)

    # Construct DataFrames
    all_loci = sorted(set(snp_map.keys()) | set(len_map.keys()) | (target_union if args.allow_extra_loci=="no" else set()))
    if not all_loci:
        sys.exit("[FATAL] No loci accumulated. Are your locus tokens consistent between targets and consensus?")
    tab_snps = pd.DataFrame(index=all_loci, columns=samples, dtype=float)
    tab_len  = pd.DataFrame(index=all_loci, columns=samples, dtype=float)

    for locus in all_loci:
        for smp in samples:
            if smp in snp_map[locus]:
                tab_snps.at[locus, smp] = snp_map[locus][smp]
            if smp in len_map[locus]:
                tab_len.at[locus, smp] = len_map[locus][smp]

    # Write outputs
    tab_snps.to_csv(out_R / "Table_SNPs.csv")
    tab_len.to_csv(out_R / "Table_consensus_length.csv")

    # Manifest: how many consensus files per sample
    man = pd.DataFrame(
        Counter([s for _, s, _ in files]).items(),
        columns=["sample", "n_consensus_files"]
    )
    # ensure all samples present
    for s in samples:
        if s not in man["sample"].values:
            man.loc[len(man)] = [s, 0]
    man.sort_values("sample").to_csv(out_R / "consensus_files_per_sample.csv", index=False)

    print(f"[OK] Wrote tables to: {out_R}")
    print(f"[INFO] loci in output: {tab_snps.shape[0]} | samples: {tab_snps.shape[1]}")

if __name__ == "__main__":
    main()
