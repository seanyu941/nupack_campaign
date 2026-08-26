"""Generate data/variants.csv.

This is a placeholder catalog so the pipeline has 400 designs to chew on
without shipping real project sequences. Replace data/variants.csv with the
actual designs and nothing downstream changes, the loader only needs the three
columns.

Sequences are drawn against a per-family GC target with two constraints that
real designs usually carry anyway: no homopolymer run longer than four, and no
GG-rich stretch that would push a G-quadruplex. Output is deterministic given
the seed.

    python scripts/make_variant_catalog.py --out data/variants.csv --seed 20260701
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np

# family name -> (count, length range, GC target)
FAMILIES: dict[str, tuple[int, tuple[int, int], float]] = {
    "toehold": (110, (28, 34), 0.48),
    "stem_loop": (110, (30, 38), 0.55),
    "linear": (90, (26, 32), 0.42),
    "branch": (90, (32, 40), 0.52),
}

HOMOPOLYMER = re.compile(r"(A{5,}|C{5,}|G{5,}|T{5,})")
G_RUN = re.compile(r"(GGG.{1,7}){3}GGG")


def passes_constraints(sequence: str) -> bool:
    return HOMOPOLYMER.search(sequence) is None and G_RUN.search(sequence) is None


def draw_sequence(rng: np.random.Generator, length: int, gc_target: float) -> str:
    """Sample until the sequence clears the constraints.

    The rejection loop is fine here because the constraints reject only a few
    percent of draws at these lengths.
    """
    weights = [
        (1 - gc_target) / 2,  # A
        gc_target / 2,        # C
        gc_target / 2,        # G
        (1 - gc_target) / 2,  # T
    ]
    for _ in range(200):
        sequence = "".join(rng.choice(list("ACGT"), size=length, p=weights))
        if passes_constraints(sequence):
            return sequence
    raise RuntimeError(
        f"could not draw a {length} nt sequence at GC {gc_target} inside 200 tries"
    )


def build_rows(seed: int) -> list[dict[str, str]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, str]] = []

    for family, (count, (min_len, max_len), gc_target) in FAMILIES.items():
        for index in range(1, count + 1):
            length = int(rng.integers(min_len, max_len + 1))
            # Spread the GC target a little so the catalog is not one narrow band.
            gc = float(np.clip(rng.normal(gc_target, 0.05), 0.30, 0.70))
            rows.append(
                {
                    "name": f"{family}_{index:03d}",
                    "sequence": draw_sequence(rng, length, gc),
                    "design_family": family,
                }
            )

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/variants.csv")
    parser.add_argument("--seed", type=int, default=20260701)
    args = parser.parse_args()

    rows = build_rows(args.seed)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        # LF rather than the csv module default of CRLF, so git diffs stay clean.
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "design_family"], lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} variants to {out}")


if __name__ == "__main__":
    main()
