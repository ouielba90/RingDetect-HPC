# RingDetect-HPC

[![PyPI version](https://badge.fury.io/py/ringdetect-hpc.svg?v=1.3.0)](https://pypi.org/project/ringdetect-hpc/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A blisteringly fast command-line utility and Python library for detecting specific molecular rings (cycles) from static structural coordinate files (XYZ, PDB, MOL, CSV), Python objects (ASE, RDKit), and massive Molecular Dynamics (MD) trajectories.

Built for computational chemistry pipelines (like EDDB, topological analysis, and materials science), RingDetect-HPC pairs a C/Python frontend for dynamic memory management with a Fortran 2008 backend for heavy graph traversal. It is deeply parallelized using OpenMP to process massive macromolecular systems and crystal lattices in fractions of a second.

## 🚀 Key Features

  * **MD Trajectory Analysis:** Natively parse multi-frame `.xyz` files (10,000+ frames) or read directly from standard input pipelines to analyze ring dynamics over time.
  * **Dual Interface:** Run it as a standalone, zero-dependency command-line executable, or `import` it natively into your Python workflows with zero-copy array pointers.
  * **Universal Periodic Table:** Automatically calculates dynamic bond cutoffs based on comprehensive covalent radii for the entire periodic table (Elements 1-96).
  * **Fragment Masking (Subsetting):** Instantly isolate specific active sites or ligands in massive proteins by passing a binary mask, entirely bypassing the O(N^2) graph math for ignored atoms.
  * **Periodic Boundary Conditions (PBC):** Seamlessly handles MOFs, zeolites, and crystal lattices by applying the Minimum Image Convention (MIC) during spatial hashing.
  * **Geometric Planarity Checking:** Automatically computes normal vectors on the fly to detect if a cycle is geometrically planar (aromaticity/conformation indicator).
  * **O(N) Spatial Hashing:** Instantly builds adjacency matrices without distance bottlenecks.
  * **Zero-Allocation Search:** Uses in-place array mutation during Depth-First Search (DFS) to completely eliminate RAM allocation locks.

-----

## 📦 Installation

You can install RingDetect-HPC in two ways, depending on your workflow.

### Option 1: Python Library (Recommended)

Pre-compiled binaries (wheels) are available for Linux and macOS via PyPI. This requires zero C/Fortran compilers on your end.

```bash
pip install ringdetect-hpc
```

*(Note: While the PyPI package is `ringdetect-hpc`, the Python module is imported as `ringdetect`.)*

### Option 2: Bare-Metal CLI (For HPC Clusters)

If you want the standalone executable, you can compile from source. **Zero dependencies required** other than standard compilers (`gcc`, `gfortran`, `make`).

```bash
git clone [https://github.com/ouielba90/RingDetect-HPC.git](https://github.com/ouielba90/RingDetect-HPC.git)
cd RingDetect-HPC
make
```

*(This generates a highly optimized executable named `ring_detector`.)*

-----

## 🐍 Python Usage

The Python API uses `ctypes` to pass data directly into Fortran's RAM without writing temporary files.

### 1. Raw Coordinates & Fragment Masking

```python
from ringdetect import find_rings

x = [0.0000, 1.2098, 1.2098, 0.0000, -1.2098, -1.2098]
y = [1.3970, 0.6985, -0.6985, -1.3970, -0.6985, 0.6985]
z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
elements = ["C", "C", "C", "C", "C", "C"]

# Find up to 6-membered rings
rings = find_rings(x, y, z, elements, max_ring=6)

# Fragment Masking: Only search for rings within a specific subset of atoms
active_subset = [0, 1, 2, 3, 4, 5]
masked_rings = find_rings(x, y, z, elements, active_atoms=active_subset)

print(rings)
# Output: [{'size': 6, 'planar': True, 'indices': [1, 2, 3, 4, 5, 6]}]
```

### 2. Third-Party Integrations (ASE & RDKit)

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

### 3. Molecular Dynamics (MDAnalysis)

Loop over massive compressed trajectories natively in Python without writing intermediate files:

```python
import MDAnalysis as mda
from ringdetect import find_rings

u = mda.Universe("protein.pdb", "trajectory.xtc")
protein = u.select_atoms("protein")

for ts in u.trajectory:
    pos = protein.positions
    rings = find_rings(pos[:,0], pos[:,1], pos[:,2], list(protein.elements))
    print(f"Frame {ts.frame}: {len(rings)} rings detected.")
```

-----

## 💻 CLI Usage

```bash
./ring_detector <molecule_file> [OPTIONS]
```
*(Note: To read from standard input, use `-` as the `<molecule_file>` name).*

### Command Line Options

| Flag | Description | Default | Example |
|------|-------------|---------|---------|
| `-h`, `--help` | Show the help menu. | | `./ring_detector -h` |
| `-f` | Set the input coordinate format (`xyz`, `raw`, `csv`, `idx`, `pdb`, `mol`). | `xyz` | `-f pdb` |
| `-a` | Active atoms mask (comma-separated ranges or single indices) to isolate sub-graphs. | `All` | `-a 1-15,30,45-50` |
| `-c` | Set unit cell dimensions (`X Y Z`) for Periodic Boundary Conditions. | `0.0 0.0 0.0`| `-c 15.5 15.5 15.5` |
| `-m` | Maximum ring depth to search (Safely capped at 100). | `6` | `-m 10` |
| `-r` | Target specific ring sizes (comma-separated). Overrides `-m`. | `All` | `-r 5,6` |
| `-p` | Number of OpenMP threads to use. (`0` = All physical cores) | `0` | `-p 8` |
| `-s` | Character used to separate atom indices in the output file. | `' '` | `-s -` |
| `-j` | Output results in strict JSON format instead of standard text. | `OFF` | `-j` |

### Example Runs

Analyze a 10,000-frame MD trajectory using 8 P-cores and outputting to JSON:
```bash
./ring_detector md_trajectory.xyz -f xyz -j -p 8
```

Pipe directly from GROMACS to avoid saving a massive text file to your disk:
```bash
gmx trjconv -s prod.tpr -f prod.xtc -o .xyz | ./ring_detector - -f xyz -j
```

-----

## 📄 Output Formats

The engine dynamically generates output using the base name of your input. For multi-frame trajectories, output is automatically grouped by frame.

### Standard Text (`.rings`)

Heavily optimized for downstream text parsing. Includes geometric planarity detection.

```text
=== FRAME 1 ===
5-MR: 1-5-4-10-9
6-MR (PLANAR): 1-2-15-14-13-12
=== FRAME 2 ===
5-MR: 1-5-4-10-9
```

### JSON Format (`-j`)

Strict JSON payload, perfect for loading directly into Pandas DataFrames for plotting frame-by-frame network stability.

```json
{
  "molecule": "md_trajectory",
  "frames": [
    {
      "frame": 1,
      "total_atoms": 60,
      "rings": [
        {"size": 5, "planar": false, "indices": [1, 5, 4, 10, 9]},
        {"size": 6, "planar": true, "indices": [1, 2, 15, 14, 13, 12]}
      ]
    }
  ]
}
```

## 🏗️ Architecture Stack

1.  **`src/main.c` / `ringdetect/engine.py`:** Executes data parsing, manages multi-frame I/O streams, translates atomic symbols to covalent radii, allocates contiguous memory blocks, applies masking logic, and passes zero-copy pointers to Fortran.
2.  **`src/ring_engine.f90`:** Receives pointers via `iso_c_binding`, builds the OpenMP spatial hash (applying Minimum Image Convention if PBC is active), and launches an optimized, lock-free Depth-First Search with in-flight vector math to isolate and classify the Minimum Cycle Basis. Temp files are streamed dynamically and merged instantly upon completion.
```
