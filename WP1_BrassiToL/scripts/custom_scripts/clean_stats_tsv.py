#This Python script helps clean the HybPiper stats files as follows.
#1) Only the first occurrence of the header row is kept;
#2) Only the final occurrence of each sample (in column "Name") is kept, removing any previous results.
#This script was created with the help of ChatGPT on 6 May 2025.
#See the following link for details: https://chatgpt.com/share/6819c037-31a4-800c-8be1-4d0129ef7ad8


#!/usr/bin/env python3

import argparse
import gzip
import os

def open_by_extension(filepath, mode="rt"):
    """Automatically choose open function based on file extension."""
    return gzip.open(filepath, mode) if filepath.endswith(".gz") else open(filepath, mode)

def clean_stats_file(stats_input_file, stats_output_file):
    last_entries = {}
    header = None

    with open_by_extension(stats_input_file, "rt") as infile:
        for line in infile:
            if line.startswith("Name"):
                if header is None:
                    header = line.strip()  # Save only the first header
                continue  # Skip repeated headers
            parts = line.strip().split("\t")
            if not parts:
                continue
            name = parts[0]
            last_entries[name] = line.strip()

    write_mode = "wt"
    with open_by_extension(stats_output_file, write_mode) as outfile:
        outfile.write(header + "\n")
        for entry in last_entries.values():
            outfile.write(entry + "\n")

def main():
    parser = argparse.ArgumentParser(description="Clean stats TSV file: remove repeated headers and keep only the last entry per Name.")
    parser.add_argument("--stats_input_file", required=True, help="Path to the input stats TSV or TSV.GZ file")
    parser.add_argument("--stats_output_file", required=True, help="Path to the cleaned output TSV or TSV.GZ file")

    args = parser.parse_args()
    clean_stats_file(args.stats_input_file, args.stats_output_file)

if __name__ == "__main__":
    main()
