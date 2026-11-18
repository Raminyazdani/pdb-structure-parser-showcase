"""Read atom locations from a single PDB model."""

import argparse
from collections import Counter
from dataclasses import dataclass
from math import dist, isfinite
from pathlib import Path
from typing import Iterable


class StructureError(ValueError):
    """The input cannot be summarized as a single set of atom locations."""


@dataclass(frozen=True)
class Atom:
    kind: str
    serial: int
    xyz: tuple[float, float, float]


def read_atoms(lines: Iterable[str]) -> dict[int, Atom]:
    atoms = {}
    seen_model = False
    for line_number, line in enumerate(lines, 1):
        line = line.rstrip("\r\n")
        kind = line[:6].strip()
        if kind == "MODEL":
            if seen_model or atoms:
                raise StructureError(f"line {line_number}: expected a single model")
            seen_model = True
        if kind not in ("ATOM", "HETATM"):
            continue
        if len(line) < 54:
            raise StructureError(f"line {line_number}: coordinates must extend through column 54")
        try:
            serial = int(line[6:11])
        except ValueError:
            raise StructureError(f"line {line_number}: invalid decimal atom serial") from None
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            raise StructureError(f"line {line_number}: invalid coordinate") from None
        if not all(isfinite(value) for value in xyz):
            raise StructureError(f"line {line_number}: non-finite coordinate")
        if serial in atoms:
            raise StructureError(f"line {line_number}: duplicate atom serial {serial}")
        atoms[serial] = Atom(kind, serial, xyz)
    if not atoms:
        raise StructureError("no ATOM or HETATM coordinates found")
    return atoms


def summarize(atoms: dict[int, Atom], pair: list[int] | None) -> str:
    counts = Counter(atom.kind for atom in atoms.values())
    summary = f"Records: {len(atoms)} (ATOM {counts['ATOM']}, HETATM {counts['HETATM']})"
    if pair is None:
        return summary
    for serial in pair:
        if serial not in atoms:
            raise StructureError(f"atom serial {serial} was not found")
    first, second = (atoms[serial] for serial in pair)
    distance = dist(first.xyz, second.xyz)
    return f"{summary}\nDistance {first.serial} -> {second.serial}: {distance:.3f} angstrom"


def main() -> None:
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("file", type=Path, help="single-model ASCII PDB file")
    cli.add_argument("--distance", metavar=("FROM", "TO"), nargs=2, type=int,
                     help="measure the separation between two atom serial numbers")
    args = cli.parse_args()
    try:
        with args.file.open(encoding="ascii") as stream:
            result = summarize(read_atoms(stream), args.distance)
    except (OSError, UnicodeError, StructureError) as error:
        cli.exit(2, f"{args.file}: {error}\n")
    print(result)


if __name__ == "__main__":
    main()
