# PDB coordinate reader

Count atom records in a single-model Protein Data Bank (PDB) coordinate file and measure the straight-line distance between two atom serial numbers. Python 3.10 or newer supplies all dependencies.

## Try the example

From the repository root, run

```sh
python3 pdb_points.py examples/points.pdb --distance 17 83
```

```text
Records: 3 (ATOM 2, HETATM 1)
Distance 17 -> 83: 13.000 angstrom
```

The example contains two fictional sites of type `ATOM` and one of type `HETATM`. Site 83 lies 3 Å along x, 4 Å along y and 12 Å along z from site 17. The squared distance is 3² + 4² + 12² = 169 Å², giving a separation of 13 Å.

## Reading another structure

Pass the path to an ASCII PDB file. Omit `--distance` for counts alone. Each coordinate record must have a decimal serial in columns 7–11 and finite x, y and z coordinates in columns 31–38, 39–46 and 47–54. Adjacent negative coordinates are read from their column positions. The reader counts each coordinate record, including each alternate location. It requires unique atom serials and accepts a single model. An input error exits with code 2. A malformed coordinate record includes its line number in the diagnostic.

## Run the tests

```sh
python3 -m unittest discover -s tests -v
```
