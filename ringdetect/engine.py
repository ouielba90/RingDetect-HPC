import ctypes
import os
import tempfile

# Locate and load the shared Fortran/C library
_lib_path = os.path.join(os.path.dirname(__file__), 'libringengine.so')
_lib = ctypes.CDLL(_lib_path)

# ---------------------------------------------------------
# Define the C-types signature to match ring_engine.f90
# ---------------------------------------------------------
_lib.find_rings.argtypes = [
    ctypes.POINTER(ctypes.c_int),       # n_atoms
    ctypes.POINTER(ctypes.c_double),    # x
    ctypes.POINTER(ctypes.c_double),    # y
    ctypes.POINTER(ctypes.c_double),    # z
    ctypes.POINTER(ctypes.c_double),    # radii
    ctypes.POINTER(ctypes.c_int),       # max_ring
    ctypes.c_char,                      # sep (PASSED BY VALUE)
    ctypes.POINTER(ctypes.c_int),       # threads
    ctypes.POINTER(ctypes.c_int),       # target_rings
    ctypes.c_char_p,                    # out_filename
    ctypes.POINTER(ctypes.c_double)     # cell (NEW: Periodic Boundaries)
]
_lib.find_rings.restype = None

def _get_covalent_radius(element):
    """Helper to map element symbols to covalent radii (in Angstroms)."""
    radii = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "PD": 1.39}
    return radii.get(element.upper(), 1.00)

def find_rings(x, y, z, elements, max_ring=6, cell=None):
    """
    Python API for the RingDetect HPC Engine.
    """
    n_atoms = len(x)
    if not (n_atoms == len(y) == len(z) == len(elements)):
        raise ValueError("Coordinates and elements lists must have the same length.")
    
    # 1. Convert Python lists to C-compatible arrays
    c_n_atoms = ctypes.c_int(n_atoms)
    c_x = (ctypes.c_double * n_atoms)(*x)
    c_y = (ctypes.c_double * n_atoms)(*y)
    c_z = (ctypes.c_double * n_atoms)(*z)
    
    radii = [_get_covalent_radius(e) for e in elements]
    c_radii = (ctypes.c_double * n_atoms)(*radii)
    
    c_max_ring = ctypes.c_int(max_ring)
    c_sep = ctypes.c_char(b',')  # Critical: Pass as byte char
    c_threads = ctypes.c_int(0)  # 0 = Use max OpenMP threads
    
    # Target all rings by default
    c_targets = (ctypes.c_int * 100)(*([1] * 100))
    
    # Handle Periodic Boundary Conditions
    if cell is None:
        cell = [0.0, 0.0, 0.0]
    c_cell = (ctypes.c_double * 3)(*cell)
    
    # Create a safe temporary file for the Fortran engine to write to
    fd, temp_path = tempfile.mkstemp(suffix=".tmp")
    os.close(fd)
    c_filename = ctypes.c_char_p(temp_path.encode('utf-8'))
    
    # 2. Execute the Fortran Engine
    _lib.find_rings(
        ctypes.byref(c_n_atoms), c_x, c_y, c_z, c_radii, 
        ctypes.byref(c_max_ring), c_sep, ctypes.byref(c_threads), 
        c_targets, c_filename, c_cell
    )
    
    # 3. Parse the results back into a Python dictionary
    rings = []
    with open(temp_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            parts = line.split(': ')
            if len(parts) != 2: continue
            
            header, indices_str = parts
            size = int(header.split('-')[0])
            is_planar = "PLANAR" in header
            indices = [int(idx) for idx in indices_str.split(',')]
            
            rings.append({
                "size": size,
                "planar": is_planar,
                "indices": indices
            })
            
    # Clean up
    os.remove(temp_path)
    return rings

def find_rings_ase(atoms, max_ring=6):
    """
    Convenience wrapper for the Atomic Simulation Environment (ASE).
    Extracts coordinates, elements, and cell dimensions directly from an ASE Atoms object.

    Parameters:
    -----------
    atoms : ase.Atoms
        The ASE Atoms object to analyze.
    max_ring : int, optional
        Maximum depth to search for rings (default is 6).

    Returns:
    --------
    list of dict
        A list of detected rings, e.g., [{'size': 6, 'planar': True, 'indices': [0,1,2,3,4,5]}]
    """
    # Extract native Python lists for the ctypes interface
    positions = atoms.get_positions()
    x = positions[:, 0].tolist()
    y = positions[:, 1].tolist()
    z = positions[:, 2].tolist()
    elements = atoms.get_chemical_symbols()

    # Extract Periodic Boundary Conditions (PBC) if they exist
    cell = [0.0, 0.0, 0.0]
    if any(atoms.pbc):
        cell_lengths = atoms.cell.lengths()
        cell = [
            float(cell_lengths[0]) if atoms.pbc[0] else 0.0,
            float(cell_lengths[1]) if atoms.pbc[1] else 0.0,
            float(cell_lengths[2]) if atoms.pbc[2] else 0.0
        ]

    return find_rings(x, y, z, elements, max_ring=max_ring, cell=cell)

def find_rings_rdkit(mol, max_ring=6, conf_id=-1):
    """
    Convenience wrapper for RDKit.
    Extracts 3D coordinates and elements directly from an RDKit Mol object.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object to analyze. Must have 3D coordinates.
    max_ring : int, optional
        Maximum depth to search for rings (default is 6).
    conf_id : int, optional
        The conformer ID to extract coordinates from (default is -1, the active conformer).

    Returns:
    --------
    list of dict
        A list of detected rings, e.g., [{'size': 6, 'planar': True, 'indices': [0,1,2,3,4,5]}]
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("RDKit molecule has no conformers! Please generate 3D coordinates (e.g., using AllChem.EmbedMolecule) before ring detection.")

    conf = mol.GetConformer(conf_id)
    n_atoms = mol.GetNumAtoms()

    x, y, z, elements = [], [], [], []

    # Extract coordinates and symbols
    for i in range(n_atoms):
        pos = conf.GetAtomPosition(i)
        x.append(pos.x)
        y.append(pos.y)
        z.append(pos.z)
        elements.append(mol.GetAtomWithIdx(i).GetSymbol())

    # RDKit standard molecules do not handle Periodic Boundary Conditions natively
    # like ASE does, so we explicitly disable PBCs by passing a zeroed cell.
    cell = [0.0, 0.0, 0.0]

    return find_rings(x, y, z, elements, max_ring=max_ring, cell=cell)
