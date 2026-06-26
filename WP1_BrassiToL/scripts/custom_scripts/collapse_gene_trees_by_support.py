#!/usr/bin/env python3
"""
Collapse internal branches in IQ-TREE gene trees based on thresholds.

- Supports node labels like:
    - "95/100"   (SH-aLRT / UFBoot)
    - "100"      (single support value, treated as UFBoot unless --prefer-alrt)
    - NHX / BEAST-style comments are ignored unless numeric is directly parseable

- Collapsing rule:
    Collapse an internal branch if ANY of these is true (default):
      (a) SH-aLRT < --alrt-threshold (if available)
      (b) UFBoot  < --ufboot-threshold (if available)
      (c) branch_length < --bl-threshold (if provided)

  Use --logic all to require that ALL available criteria fail before collapsing.

Outputs one collapsed tree per input, preserving leaf names and writing Newick.

Author: Kasper’s helper
"""

import argparse, glob, os, re, sys
from typing import Optional, Tuple

try:
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree, Clade
except Exception as e:
    sys.stderr.write("[ERR] Biopython not available. Please install with: pip install biopython\n")
    raise

PAIR_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)[/\|](\d+(?:\.\d+)?)\s*$')  # e.g., "93/100" or "93|100"
NUM_RE  = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*$')                      # e.g., "100" or "95.5"

def parse_support(label: Optional[str], prefer_alrt: bool=False) -> Tuple[Optional[float], Optional[float]]:
    """
    Try to parse 'aLRT/UFBoot' or a single number out of a node label.
    Returns (alrt, ufboot) in percentages if present, else None.
    """
    if label is None:
        return (None, None)

    s = str(label)
    m = PAIR_RE.match(s)
    if m:
        a, b = m.group(1), m.group(2)
        try:
            return (float(a), float(b))
        except:
            return (None, None)

    m = NUM_RE.match(s)
    if m:
        v = float(m.group(1))
        if prefer_alrt:
            return (v, None)
        else:
            return (None, v)

    # Try to rescue from comments like [&SH_ALRT=99,UFBOOT=100] – optional
    alrt = None
    ufbt = None
    if "SH" in s or "ALRT" in s or "UF" in s or "boot" in s.lower():
        try:
            if "ALRT" in s.upper():
                m2 = re.search(r'ALRT\s*=\s*([\d.]+)', s, flags=re.IGNORECASE)
                if m2: alrt = float(m2.group(1))
            if "UF" in s.upper() or "BOOT" in s.lower():
                m3 = re.search(r'(UFBOOT|UF|BOOTSTRAP)\s*=\s*([\d.]+)', s, flags=re.IGNORECASE)
                if m3: ufbt = float(m3.group(2))
        except:
            pass
    return (alrt, ufbt)

def should_collapse(clade: Clade, alrt_th: float, ufboot_th: float, bl_th: Optional[float],
                    logic: str, prefer_alrt: bool) -> bool:
    """
    Decide if the edge above 'clade' should be collapsed.
    We look at the label on 'clade' itself (standard Newick puts support on the child clade).
    """
    # Skip root: there is no parent edge to collapse
    # In Bio.Phylo, root often has no parent; we detect by absence of clade.branch_length at parent edge.
    # Here, we treat the clade itself; collapsing will be done on this clade in the Tree object.

    alrt, ufbt = parse_support(getattr(clade, "name", None), prefer_alrt=prefer_alrt)

    tests = []
    if alrt is not None and alrt_th is not None:
        tests.append(alrt < alrt_th)
    if ufbt is not None and ufboot_th is not None:
        tests.append(ufbt < ufboot_th)
    if bl_th is not None and clade.branch_length is not None:
        try:
            tests.append(float(clade.branch_length) < bl_th)
        except:
            pass

    if not tests:
        return False  # nothing to decide on this edge

    if logic == "all":
        return all(tests)
    else:
        return any(tests)

def main():
    ap = argparse.ArgumentParser(description="Collapse IQ-TREE gene trees by SH-aLRT/UFBoot and/or branch length thresholds.")
    ap.add_argument("--input-glob", required=True, help="Glob for input .tree/.treefile (one tree per file).")
    ap.add_argument("--outdir",     required=True, help="Output directory for collapsed trees.")
    ap.add_argument("--alrt-threshold",  type=float, default=70.0, help="SH-aLRT threshold (percent).")
    ap.add_argument("--ufboot-threshold",type=float, default=80.0, help="UFBoot threshold (percent).")
    ap.add_argument("--bl-threshold",    type=float, default=None, help="Branch length threshold; collapse if BL < this.")
    ap.add_argument("--logic", choices=["any","all"], default="any",
                    help="Collapse if ANY test fails (default) or require ALL to fail.")
    ap.add_argument("--prefer-alrt", action="store_true",
                    help="If only a single numeric label is found, interpret it as SH-aLRT (default: as UFBoot).")
    ap.add_argument("--suffix", default="_collapsed.tree", help="Suffix for output filenames.")
    args = ap.parse_args()

    files = sorted(glob.glob(args.input_glob))
    if not files:
        sys.stderr.write(f"[WARN] No files matched: {args.input_glob}\n")

    os.makedirs(args.outdir, exist_ok=True)
    n_in = 0
    n_out = 0

    for fp in files:
        n_in += 1
        try:
            tree: Tree = Phylo.read(fp, "newick")
        except Exception as e:
            sys.stderr.write(f"[WARN] Skipping unreadable tree: {fp} ({e})\n")
            continue

        # Collect clades to collapse (internal nodes only; skip terminals)
        to_collapse = []
        for cl in tree.find_clades(order="level"):
            if cl.is_terminal():
                continue
            # Root guard: collapsing root makes no sense
            # Heuristic: Bio.Phylo lacks parent links; we avoid collapsing the very first internal node if it appears to be root
            # This is safe enough because thresholds usually target low-support internal edges deeper in the tree.
            if cl == tree.root:
                continue

            if should_collapse(cl, args.alrt_threshold, args.ufboot_threshold, args.bl_threshold, args.logic, args.prefer_alrt):
                to_collapse.append(cl)

        # Collapse marked clades (from deepest to shallower to avoid invalidating references)
        # Using the Tree.collapse(clade) API.
        for cl in sorted(to_collapse, key=lambda c: c.count_terminals(), reverse=True):
            try:
                tree.collapse(cl)
            except Exception as e:
                sys.stderr.write(f"[WARN] Could not collapse a clade in {os.path.basename(fp)}: {e}\n")

        outname = os.path.basename(fp)
        if outname.lower().endswith(".treefile"):
            outname = outname[:-9] + args.suffix
        elif outname.lower().endswith(".tree"):
            outname = outname[:-5] + args.suffix
        else:
            outname = outname + args.suffix

        outpath = os.path.join(args.outdir, outname)
        try:
            Phylo.write(tree, outpath, "newick")
            n_out += 1
        except Exception as e:
            sys.stderr.write(f"[ERR] Failed to write {outpath}: {e}\n")

    sys.stderr.write(f"[OK] Collapsed {n_out}/{n_in} trees into {args.outdir}\n")

if __name__ == "__main__":
    main()
