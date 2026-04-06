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
    ctypes.POINTER(ctypes.c_double),    # cell 
    ctypes.POINTER(ctypes.c_int)        # active_mask
]
_lib.find_rings.restype = None

def _get_covalent_radius(element):
    """Helper to map element symbols to covalent radii (in Angstroms)."""
    radii = {
        "H": 0.31, "HE": 0.28, "LI": 1.28, "BE": 0.96, "B": 0.84, "C": 0.76, 
        "N": 0.71, "O": 0.66, "F": 0.57, "NE": 0.58, "NA": 1.66, "MG": 1.41, 
        "AL": 1.21, "SI": 1.11, "P": 1.07, "S": 1.05, "CL": 1.02, "AR": 1.06, 
        "K": 2.03, "CA": 1.76, "SC": 1.70, "TI": 1.60, "V": 1.53, "CR": 1.39, 
        "MN": 1.39, "FE": 1.32, "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22, 
        "GA": 1.22, "GE": 1.20, "AS": 1.19, "SE": 1.20, "BR": 1.20, "KR": 1.16, 
        "RB": 2.20, "SR": 1.95, "Y": 1.90, "ZR": 1.75, "NB": 1.64, "MO": 1.54, 
        "TC": 1.47, "RU": 1.46, "RH": 1.42, "PD": 1.39, "AG": 1.45, "CD": 1.44, 
        "IN": 1.42, "SN": 1.39, "SB": 1.39, "TE": 1.38, "I": 1.39, "XE": 1.40, 
        "CS": 2.44, "BA": 2.15, "LA": 2.07, "CE": 2.04, "PR": 2.03, "ND": 2.01, 
        "PM": 1.99, "SM": 1.98, "EU": 1.98, "GD": 1.96, "TB": 1.94, "DY": 1.92, 
        "HO": 1.92, "ER": 1.89, "TM": 1.90, "YB": 1.87, "LU": 1.87, "HF": 1.75, 
        "TA": 1.70, "W": 1.62, "RE": 1.51, "OS": 1.44, "IR": 1.41, "PT": 1.36, 
        "AU": 1.36, "HG": 1.32, "TL": 1.45, "PB": 1.46, "BI": 1.48, "PO": 1.40, 
        "AT": 1.50, "RN": 1.50, "FR": 2.60, "RA": 2.21, "AC": 2.15, "TH": 2.06, 
        "PA": 2.00, "U": 1.96, "NP": 1.90, "PU": 1.87, "AM": 1.80, "CM": 1.69
    }
    # Clean string format to uppercase for safe dictionary matching
    clean_elem = str(element).strip().upper()
    return radii.get(clean_elem, 1.00)

def find_rings(x, y, z, elements, max_ring=6, cell=None, active_atoms=None):
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
    c_sep = ctypes.c_char(b',')  
    c_threads = ctypes.c_int(0)  
    c_targets = (ctypes.c_int * 100)(*([1] * 100))
    
    if cell is None:
        cell = [0.0, 0.0, 0.0]
    c_cell = (ctypes.c_double * 3)(*cell)
    
    # --- NEW: Handle Active Atoms (Fragment Masking) ---
    if active_atoms is None:
        mask = [1] * n_atoms
    else:
        mask = [0] * n_atoms
        for idx in active_atoms:
            # Python uses 0-based indexing!
            if 0 <= idx < n_atoms:
                mask[idx] = 1
    c_active_mask = (ctypes.c_int * n_atoms)(*mask)
    # ---------------------------------------------------

    fd, temp_path = tempfile.mkstemp(suffix=".tmp")
    os.close(fd)
    c_filename = ctypes.c_char_p(temp_path.encode('utf-8'))
    
    # 2. Execute the Fortran Engine (Pass c_active_mask!)
    _lib.find_rings(
        ctypes.byref(c_n_atoms), c_x, c_y, c_z, c_radii, 
        ctypes.byref(c_max_ring), c_sep, ctypes.byref(c_threads), 
        c_targets, c_filename, c_cell, c_active_mask
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
            
    os.remove(temp_path)
    return rings

def find_rings_ase(atoms, max_ring=6, active_atoms=None):
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

def find_rings_rdkit(mol, max_ring=6, conf_id=-1, active_atoms=None):
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

    return find_rings(x, y, z, elements, max_ring=max_ring, cell=cell, active_atoms=active_atoms)
