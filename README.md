# RingDetect-HPC

[![PyPI version](https://badge.fury.io/py/ringdetect.svg)](https://pypi.org/project/ringdetect/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A blisteringly fast command-line utility and Python library for detecting specific molecular rings (cycles) from various structural coordinate files (XYZ, PDB, MOL, CSV) and Python objects (ASE, RDKit). 

Built for computational chemistry pipelines (like EDDB, topological analysis, and materials science), RingDetect-HPC pairs a C/Python frontend for dynamic memory management with a Fortran 2008 backend for heavy graph traversal. It is deeply parallelized using OpenMP to process massive macromolecular systems and crystal lattices in fractions of a second.

## 🚀 Key Features
* **Dual Interface:** Run it as a standalone, zero-dependency command-line executable, or `import` it natively into your Python workflows with zero-copy array pointers.
* **Periodic Boundary Conditions (PBC):** Seamlessly handles MOFs, zeolites, and crystal lattices by applying the Minimum Image Convention (MIC) during spatial hashing.
* **Geometric Planarity Checking:** Automatically computes normal vectors on the fly to detect if a cycle is geometrically planar (aromaticity/conformation indicator).
* **O(N) Spatial Hashing:** Instantly builds adjacency matrices without O(N^2) distance bottlenecks.
* **Zero-Allocation Search:** Uses in-place array mutation during Depth-First Search (DFS) to completely eliminate RAM allocation locks.
* **JSON Export:** Output raw paths to standard text files, or strict JSON formats for easy parsing by Pandas and web dashboards.

---

## 📦 Installation

You can install RingDetect-HPC in two ways, depending on your workflow.

### Option 1: Python Library (Recommended)
Pre-compiled binaries (wheels) are available for Linux and macOS. This requires zero C/Fortran compilers on your end.

```bash
pip install ringdetect
```

### Option 2: Bare-Metal CLI (For HPC Clusters)
If you want the standalone executable, you can compile from source. **Zero dependencies required** other than standard compilers (`gcc`, `gfortran`, `make`).

```bash
git clone [https://github.com/ouielba90/RingDetect-HPC.git](https://github.com/ouielba90/RingDetect-HPC.git)
cd RingDetect-HPC
make
```
*(This generates a highly optimized executable named `ring_detector`.)*

---

## 🐍 Python Usage

The Python API uses `ctypes` to pass data directly into Fortran's RAM without writing temporary files.

### 1. Raw Coordinates
```python
from ringdetect import find_rings

x = [0.0000, 1.2098, 1.2098, 0.0000, -1.2098, -1.2098]
y = [1.3970, 0.6985, -0.6985, -1.3970, -0.6985, 0.6985]
z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
elements = ["C", "C", "C", "C", "C", "C"]

# Find up to 6-membered rings using 4 CPU cores
rings = find_rings(x, y, z, elements, max_ring=6, threads=4)

print(rings)
# Output: [{'size': 6, 'planar': True, 'indices': [1, 2, 3, 4, 5, 6]}]
```

### 2. Third-Party Integrations (ASE & RDKit)
If you are already using popular computational chemistry libraries, you can pass those objects directly:

```python
from ringdetect import find_rings_ase, find_rings_rdkit
from ase.io import read
from rdkit import Chem

# With ASE (Automatically detects Periodic Boundary Conditions!)
atoms = read("crystal_lattice.cif")
rings = find_rings_ase(atoms, max_ring=8)

# With RDKit
mol = Chem.MolFromMolFile("molecule.mol")
rings = find_rings_rdkit(mol, max_ring=6)
```

---

## 💻 CLI Usage

```bash
./ring_detector <molecule_file> [OPTIONS]
```

### Command Line Options

| Flag | Description | Default | Example |
|------|-------------|---------|---------|
| `-h`, `--help` | Show the help menu. | | `./ring_detector -h` |
| `-f` | Set the input coordinate format (`xyz`, `raw`, `csv`, `idx`, `pdb`, `mol`). | `xyz` | `-f pdb` |
| `-c` | Set unit cell dimensions (`X Y Z`) for Periodic Boundary Conditions. | `0.0 0.0 0.0`| `-c 15.5 15.5 15.5` |
| `-m` | Maximum ring depth to search. | `6` | `-m 10` |
| `-r` | Target specific ring sizes (comma-separated). Overrides `-m`. | `All` | `-r 5,6` |
| `-p` | Number of OpenMP threads to use. (`0` = All physical cores) | `0` | `-p 8` |
| `-s` | Character used to separate atom indices in the output file. | `' '` | `-s -` |
| `-j` | Output results in strict JSON format instead of standard text. | `OFF` | `-j` |

### Example Run

To search a crystal XYZ file for only 5-membered and 6-membered rings, applying a 15.0 Å unit cell, using 8 CPU cores, and outputting to JSON:

```bash
./ring_detector crystal.xyz -c 15.0 15.0 15.0 -r 5,6 -p 8 -j
```

---

## 📄 Output Formats

The engine dynamically generates output using the base name of your input (e.g., running `c60.xyz` generates `c60.rings` or `c60.json` in the execution directory).

### Standard Text (`.rings`)
Heavily optimized for downstream text parsing. Includes geometric planarity detection.
```text
5-MR: 1-5-4-10-9
6-MR (PLANAR): 1-2-15-14-13-12
6-MR: 2-3-17-16-15-1
```

### JSON Format (`-j`)
Strict JSON payload, perfect for loading directly into Pandas DataFrames or web interfaces.
```json
{
  "molecule": "c60",
  "total_atoms": 60,
  "rings": [
    {"size": 5, "planar": false, "indices": [1, 5, 4, 10, 9]},
    {"size": 6, "planar": true, "indices": [1, 2, 15, 14, 13, 12]}
  ]
}
```

## 🏗️ Architecture Stack

1.  **`src/main.c` / `ringdetect/engine.py`:** Executes data parsing (skipping PDB headers, extracting ASE geometry, etc.), translates atomic symbols to covalent radii, allocates contiguous memory blocks, and passes zero-copy pointers to Fortran.
2.  **`src/ring_engine.f90`:** Receives pointers via `iso_c_binding`, builds the OpenMP spatial hash (applying Minimum Image Convention if PBC is active), and launches an optimized, lock-free Depth-First Search with in-flight vector math to isolate and classify the Minimum Cycle Basis. Temp files are streamed dynamically and merged instantly upon completion.
```
