# RingDetect-HPC

[![PyPI version](https://badge.fury.io/py/ringdetect-hpc.svg?v=1.4.0)](https://pypi.org/project/ringdetect-hpc/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A blisteringly fast command-line utility and Python library for detecting specific molecular rings (cycles) from static structural coordinate files (XYZ, PDB, MOL, CSV), Python objects (ASE, RDKit), and massive Molecular Dynamics (MD) trajectories.

Built for computational chemistry pipelines (like EDDB, topological analysis, and materials science), RingDetect-HPC pairs a C/Python frontend for dynamic memory management with a Fortran 2008 backend for heavy graph traversal. It is deeply parallelized using OpenMP to process massive macromolecular systems and crystal lattices in fractions of a second.

## 🚀 Key Features

  * **Big Data Streaming (NEW):** Generate flat CSV streams capable of processing 80,000+ frames in under 40 seconds, perfectly optimized for Apache Parquet compression.
  * **MD Trajectory Analysis:** Natively parse multi-frame `.xyz` files or read directly from standard input pipelines to analyze ring dynamics over time.
  * **Dual Interface:** Run it as a standalone, zero-dependency command-line executable, or `import` it natively into your Python workflows with zero-copy array pointers.
  * **Universal Periodic Table:** Automatically calculates dynamic bond cutoffs based on comprehensive covalent radii for the entire periodic table (Elements 1-96).
  * **Fragment Masking (Subsetting):** Instantly isolate specific active sites or ligands in massive proteins by passing a binary mask, entirely bypassing the O(N^2) graph math for ignored atoms.
  * **Periodic Boundary Conditions (PBC):** Seamlessly handles MOFs, zeolites, and crystal lattices by applying the Minimum Image Convention (MIC).
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
| `-a` | Active atoms mask (comma-separated ranges or single indices). | `All` | `-a 1-15,30,45-50` |
| `-c` | Set unit cell dimensions (`X Y Z`) for Periodic Boundary Conditions. | `0.0 0.0 0.0`| `-c 15.5 15.5 15.5` |
| `-m` | Maximum ring depth to search (Safely capped at 100). | `6` | `-m 10` |
| `-r` | Target specific ring sizes (comma-separated). Overrides `-m`. | `All` | `-r 5,6` |
| `-p` | Number of OpenMP threads to use. (`0` = All physical cores) | `0` | `-p 8` |
| `-j` | Output results in strict JSON format. | `OFF` | `-j` |
| `-v` | Output results in flat CSV format (Big Data optimized). | `OFF` | `-v` |

-----

## 📄 Output Formats & Big Data Workflow

The engine dynamically generates output using the base name of your input. 

### 1. Flat CSV (`-v`) & Apache Parquet (Big Data)

For massive multi-frame MD trajectories, JSON becomes a bottleneck. The `-v` flag tells the C-engine to stream a flat CSV at maximum speed (~2,000 frames/sec on modern hardware).

```csv
Frame,RingSize,Planar,Indices
1,6,True,1-2-15-14-13-12
1,5,False,1-5-4-10-9
```

You can then instantly compress this into an **Apache Parquet** file using `polars` in Python (yielding ~94% file size reduction and microsecond load times):

```python
import polars as pl
# See /examples/parquet/ for the full script!
pl.read_csv("trajectory.csv").write_parquet("trajectory.parquet")
```

### 2. JSON Format (`-j`)

Strict JSON payload, perfect for loading directly into web interfaces or smaller Pandas DataFrames.

```json
{
  "molecule": "md_trajectory",
  "frames": [
    {
      "frame": 1,
      "total_atoms": 60,
      "rings": [
        {"size": 5, "planar": false, "indices": [1, 5, 4, 10, 9]}
      ]
    }
  ]
}
```

### 3. Standard Text (`.rings`)

The default output. Heavily optimized for downstream text parsing.

```text
=== FRAME 1 ===
5-MR: 1-5-4-10-9
6-MR (PLANAR): 1-2-15-14-13-12
```

## 🏗️ Architecture Stack

1.  **`src/main.c` / `ringdetect/engine.py`:** Executes data parsing, manages multi-frame I/O streams, translates atomic symbols to covalent radii, allocates contiguous memory blocks, applies masking logic, and passes zero-copy pointers to Fortran.
2.  **`src/ring_engine.f90`:** Receives pointers via `iso_c_binding`, builds the OpenMP spatial hash (applying Minimum Image Convention if PBC is active), and launches an optimized, lock-free Depth-First Search with in-flight vector math to isolate and classify the Minimum Cycle Basis. Temp files are streamed dynamically and merged instantly upon completion.
```

---

### 📝 Release Notes for GitHub (v1.4.0)

**Title:** Release v1.4.0: The Big Data Update (CSV Streaming & Parquet Support) 📊

**Description:**
This release addresses the primary I/O bottleneck encountered when analyzing massive Molecular Dynamics trajectories (80,000+ frames). We have introduced a high-speed CSV streaming pipeline tailored specifically for modern Big Data architectures like Apache Parquet.

**✨ New Features:**
* **Flat CSV Streaming (`-v` flag):** Added the `-v` (values) flag to bypass JSON formatting overhead. The C-engine now streams a perfectly flat `Frame,RingSize,Planar,Indices` table directly to disk.
* **Blistering Trajectory Speed:** With the new CSV pipeline, the OpenMP engine can process and write **~2,000 frames per second** (tested on an Intel Core i7-14700HX for a 4,000-atom system).
* **Parquet Workflow Integration:** Added the `examples/parquet/` directory to demonstrate how to achieve **94% file size compression** and microsecond load times by piping the CSV output into Python's `polars` library.

**🛠️ Upgrades & Fixes:**
* Mutually exclusive output guards in the CLI parser to prevent JSON (`-j`) and CSV (`-v`) from clashing.
* Hardened internal string processing to ensure indices containing accidental commas in memory safely fallback to dashes in CSV mode.

