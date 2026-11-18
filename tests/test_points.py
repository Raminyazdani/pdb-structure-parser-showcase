import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from pdb_points import StructureError, read_atoms

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "points.pdb"


def invoke(path, *options):
    return subprocess.run(
        [sys.executable, str(ROOT / "pdb_points.py"), str(path), *options],
        cwd=path.parent, capture_output=True, text=True,
    )


class PointReaderTests(unittest.TestCase):
    def test_example_records_and_coordinates(self):
        with EXAMPLE.open(encoding="ascii") as stream:
            atoms = read_atoms(stream)
        self.assertEqual([(a.serial, a.kind, a.xyz) for a in atoms.values()], [
            (17, "ATOM", (-2.125, -0.875, 1.25)),
            (29, "ATOM", (-0.5, -4.25, -1.75)),
            (83, "HETATM", (0.875, 3.125, 13.25)),
        ])

    def test_demo_against_pencil_and_paper_geometry(self):
        # Two ATOM rows and one HETATM row. The chosen vector is (3, 4, 12).
        # Its squared length is 9 + 16 + 144 = 169.
        result = invoke(EXAMPLE, "--distance", "17", "83")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), [
            "Records: 3 (ATOM 2, HETATM 1)",
            "Distance 17 -> 83: 13.000 angstrom",
        ])

    def test_truncated_z_has_a_line_number_and_no_partial_output(self):
        header, first, second, *_ = EXAMPLE.read_text(encoding="ascii").splitlines()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "partial.pdb"
            path.write_text("\n".join([header, first, second[:53]]) + "\n", encoding="ascii")
            result = invoke(path)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, "")
        self.assertIn("line 3: coordinates must extend through column 54", result.stderr)

    def test_adjacent_negative_coordinates(self):
        line = "HETATM  321".ljust(30) + "-123.875-234.125-345.625"
        self.assertEqual(read_atoms([line])[321].xyz, (-123.875, -234.125, -345.625))

    def test_invalid_coordinate(self):
        line = "ATOM      1".ljust(30) + "  broken   1.125   4.125"
        with self.assertRaisesRegex(StructureError, "line 1: invalid coordinate"):
            read_atoms([line])

    def test_second_model_rejected(self):
        with self.assertRaisesRegex(StructureError, "line 3: expected a single model"):
            read_atoms(["MODEL        1", "ENDMDL", "MODEL        2"])

    def test_unknown_serial(self):
        result = invoke(EXAMPLE, "--distance", "17", "91")
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, "")
        self.assertIn("atom serial 91 was not found", result.stderr)


if __name__ == "__main__":
    unittest.main()
