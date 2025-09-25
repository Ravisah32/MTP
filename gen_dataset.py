#!/usr/bin/env python3
import argparse
import csv
import random
from pathlib import Path
from collections import OrderedDict

DEFAULT_PALETTE = [
    "red", "blue", "green", "yellow", "black", "white", "orange", "purple",
    "brown", "pink", "gray", "cyan", "magenta", "teal", "maroon", "olive",
    "navy", "lime", "indigo", "violet"
]

def build_palette(x: int) -> list[str]:
    if x < 1:
        raise ValueError("x (number of colors) must be >= 1")
    if x <= len(DEFAULT_PALETTE):
        return DEFAULT_PALETTE[:x]
    # deterministically extend if someone asks for more colors than we ship
    extra = [f"color_{i}" for i in range(len(DEFAULT_PALETTE)+1, x+1)]
    return DEFAULT_PALETTE + extra

def header_for(colors: list[str]) -> list[str]:
    # index, color, then TitleCase + 'Count' for each color in order
    count_cols = [f"{c.capitalize()}Count" for c in colors]
    return ["index", "color"] + count_cols

def main():
    ap = argparse.ArgumentParser(
        description="Generate a dataset of colors with prefix-sum counts per color."
    )
    ap.add_argument("n", type=int, help="Number of rows (tuples) to generate.")
    ap.add_argument("x", type=int, help="Number of distinct colors to use.")
    ap.add_argument("-o", "--output", required=True, help="Output CSV filename.")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)
    if args.n < 1:
        raise SystemExit("n must be >= 1")
    if args.x < 1:
        raise SystemExit("x must be >= 1")

    colors = build_palette(args.x)
    counts = {c: 0 for c in colors}
    fieldnames = header_for(colors)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(1, args.n + 1):
            chosen = random.choice(colors)
            counts[chosen] += 1

            row = OrderedDict()
            row["index"] = i
            row["color"] = chosen
            # prefix sums for every color in palette order
            for c in colors:
                row[f"{c.capitalize()}Count"] = counts[c]

            writer.writerow(row)

    print(f"Wrote {args.n} rows with {args.x} colors to {out_path}")

if __name__ == "__main__":
    main()
