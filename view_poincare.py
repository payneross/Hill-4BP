#!/usr/bin/env python3
"""Interactive, infinitely-zoomable Poincare section viewer.

Reads the full-precision CSV (not the quantized SVG) and renders it in an
interactive matplotlib window -- the equivalent of a MATLAB .fig: the raw
double-precision points are kept in memory and re-rasterized at every zoom
level, so detail keeps resolving as you scroll in.

Usage:
    python3 view_poincare.py                      # defaults to hill4bp_capd.csv
    python3 view_poincare.py somefile.csv
    python3 view_poincare.py somefile.csv --png out.png   # also dump a high-DPI PNG
    python3 view_poincare.py somefile.csv --pdf out.pdf   # also dump a vector PDF

Controls in the window:
    - scroll wheel / magnifier tool : zoom
    - pan tool                      : drag
    - 'home' button                 : reset view
    - save button                   : export the current view at high DPI
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    args = [a for a in sys.argv[1:]]
    path = "hill4bp_capd.csv"
    png_out = None
    pdf_out = None
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--png":
            png_out = args[i + 1]; i += 2; continue
        if a == "--pdf":
            pdf_out = args[i + 1]; i += 2; continue
        path = a; i += 1

    data = np.loadtxt(path, delimiter=",")
    x, y = data[:, 0], data[:, 1]
    print(f"Loaded {len(x):,} points from {path}")

    fig, ax = plt.subplots(figsize=(10, 8))
    # s is the marker AREA in points^2; small fixed dots stay crisp as you
    # zoom because matplotlib re-draws them from the source coords each time.
    ax.scatter(x, y, s=1.2, c="blue", linewidths=0, marker=".")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(path)
    fig.tight_layout()

    if png_out:
        fig.savefig(png_out, dpi=600)
        print(f"Wrote {png_out} at 600 dpi")
    if pdf_out:
        fig.savefig(pdf_out)  # vector: infinite zoom in any PDF reader
        print(f"Wrote {pdf_out} (vector)")

    if not (png_out or pdf_out):
        plt.show()

if __name__ == "__main__":
    main()
