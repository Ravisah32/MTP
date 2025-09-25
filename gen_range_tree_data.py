#!/usr/bin/env python3
import csv
import itertools
import random
import sys

DEFAULT_COLOR_SEQUENCE = [
    "Red", "Blue", "Green", "Yellow", "Black",
    "Orange", "Purple", "Cyan", "Magenta", "Brown",
    "Pink", "Gray", "Teal", "Maroon", "Olive",
    "Navy", "Lime", "Aqua", "Silver", "Gold",
]

def main():
    if len(sys.argv) != 4:
        print("Usage: python gen_pairs.py <tuplesize> <numberofcolor> <outputfilename>")
        sys.exit(1)

    try:
        n = int(sys.argv[1])            # tuplesize
        x = int(sys.argv[2])            # numberofcolor
        out_path = sys.argv[3]          # outputfilename
    except Exception:
        print("Error: tuplesize and numberofcolor must be integers.")
        sys.exit(1)

    if n <= 0 or x <= 1:
        print("Error: tuplesize must be > 0 and numberofcolor must be > 1.")
        sys.exit(1)

    # Build color list according to your mapping:
    # 3 -> Red, Blue, Green; 4 adds Yellow; 5 adds Black; then continue.
    colors = []
    if x <= len(DEFAULT_COLOR_SEQUENCE):
        colors = DEFAULT_COLOR_SEQUENCE[:x]
    else:
        colors = DEFAULT_COLOR_SEQUENCE[:]  # start with defaults
        # Add generic names if we run out
        for i in range(len(DEFAULT_COLOR_SEQUENCE) + 1, x + 1):
            colors.append(f"Color{i}")

    # Unordered pairs in **input order** (not lexicographic):
    # For [Red, Blue, Green] → (Red,Blue), (Red,Green), (Blue,Green)
    pairs = [(colors[i], colors[j]) for i in range(len(colors)) for j in range(i + 1, len(colors))]
    pair_col_names = [f"{a}-{b}" for (a, b) in pairs]

    header = ["index", "color"] + pair_col_names

    # Running signed difference per pair, start at 0.
    pair_values = { (a, b): 0 for (a, b) in pairs }

    # Simple uniform sampler over colors; replace with weighted sampling if needed.
    def draw_color():
        return random.choice(colors)

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)

        for i in range(1, n + 1):
            c = draw_color()

            # Update signed difference for each pair:
            # if row color == A: +1; if == B: -1; else unchanged
            for (a, b) in pairs:
                if c == a:
                    pair_values[(a, b)] += 1
                elif c == b:
                    pair_values[(a, b)] -= 1

            # Emit row
            row = [i, c] + [pair_values[(a, b)] for (a, b) in pairs]
            w.writerow(row)

    print(f"Wrote {n} rows, {x} colors, {len(pairs)} pair columns → {out_path}")

if __name__ == "__main__":
    main()
