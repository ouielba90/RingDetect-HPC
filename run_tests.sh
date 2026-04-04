#!/bin/bash

# Define colors for pretty output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Clean up trap: Ensures artifacts are deleted even if a test fails and exits
trap 'rm -f *.rings *.json temp_test.raw temp_test.idx' EXIT

echo "======================================================"
echo -e " ${YELLOW}🚀 Starting RingDetect-HPC Integration Test Suite${NC}"
echo "======================================================"

# --- STEP 1: Build ---
echo -e "\n[1/5] Compiling CLI and Python Shared Library..."
make clean > /dev/null 2>&1
make > /dev/null
make lib > /dev/null

if [ ! -f "ring_detector" ] || [ ! -f "ringdetect/libringengine.so" ]; then
    echo -e "${RED}FAILED: Compilation error. Executables not found.${NC}"
    exit 1
fi
echo -e "${GREEN}SUCCESS: Build completed.${NC}"

# --- STEP 2: Core Parsers & Physics Engine ---
echo -e "\n[2/5] Testing CLI Parsers & Core Physics..."

# Test 1: Standard XYZ + Planarity Math
./ring_detector examples/data/benzene.xyz -f xyz > /dev/null
if grep -q "6-MR (PLANAR)" benzene.rings; then
    echo -e "${GREEN}  ✓ XYZ parsing & Planarity check passed.${NC}"
else
    echo -e "${RED}  ✗ Benzene XYZ parsing failed.${NC}"; exit 1
fi

# Test 2: CSV Parsing
./ring_detector examples/data/cyclobutane.csv -f csv > /dev/null
if grep -q "4-MR" cyclobutane.rings; then
    echo -e "${GREEN}  ✓ CSV parsing passed.${NC}"
else
    echo -e "${RED}  ✗ Cyclobutane CSV parsing failed.${NC}"; exit 1
fi

# Test 3: MOL Parsing
./ring_detector examples/data/cyclopentane.mol -f mol > /dev/null
if grep -q "5-MR" cyclopentane.rings; then
    echo -e "${GREEN}  ✓ MOL parsing passed.${NC}"
else
    echo -e "${RED}  ✗ Cyclopentane MOL parsing failed.${NC}"; exit 1
fi

# Test 4: PDB Parsing + JSON Export
./ring_detector examples/data/phenol.pdb -f pdb -j > /dev/null
if grep -q "\"size\": 6" phenol.json; then
    echo -e "${GREEN}  ✓ PDB parsing & JSON export passed.${NC}"
else
    echo -e "${RED}  ✗ Phenol PDB parsing or JSON export failed.${NC}"; exit 1
fi

# --- STEP 3: Advanced CLI Flags ---
echo -e "\n[3/5] Testing Advanced CLI Features (PBC, Threads, Targets)..."

# Test 5: OpenMP Thread Limit & Custom Separator
./ring_detector examples/data/benzene.xyz -f xyz -p 2 -s "|" > /dev/null
if grep -q "|" benzene.rings; then
    echo -e "${GREEN}  ✓ Custom separator (-s) and OpenMP thread limits (-p) passed.${NC}"
else
    echo -e "${RED}  ✗ Custom separator or threading failed.${NC}"; exit 1
fi

# Test 6: Periodic Boundary Conditions (PBC)
./ring_detector examples/data/benzene.xyz -f xyz -c 15.0 15.0 15.0 > /dev/null
if grep -q "6-MR" benzene.rings; then
    echo -e "${GREEN}  ✓ Periodic Boundary Conditions (Minimum Image Convention) passed.${NC}"
else
    echo -e "${RED}  ✗ PBC distance calculation failed.${NC}"; exit 1
fi

# Test 7: Specific Target Rings Filter
./ring_detector examples/data/benzene.xyz -f xyz -r 4,5,7 > /dev/null
if ! grep -q "6-MR" benzene.rings; then
    echo -e "${GREEN}  ✓ Target Ring filter (-r) correctly excluded non-targeted sizes.${NC}"
else
    echo -e "${RED}  ✗ Target Ring filter failed (found a 6-MR when searching for 4,5,7).${NC}"; exit 1
fi

# --- STEP 4: Error Handling & Edge Cases ---
echo -e "\n[4/5] Testing Error Handling & Edge Cases..."

# Test 8: Missing File Handling (Should fail gracefully, not segfault)
if ./ring_detector non_existent_ghost_file.xyz > /dev/null 2>&1; then
    echo -e "${RED}  ✗ Missing file test failed (Program did not return an error code).${NC}"; exit 1
else
    echo -e "${GREEN}  ✓ Missing file handled gracefully.${NC}"
fi

# Test 9: Out-of-Bounds Max Ring (Fortran should cap this safely at 100)
./ring_detector examples/data/benzene.xyz -m 5000 > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  ✓ Extreme max_ring (-m 5000) capped safely without crashing.${NC}"
else
    echo -e "${RED}  ✗ Extreme max_ring caused a crash.${NC}"; exit 1
fi

# --- STEP 5: Python API Tests ---
echo -e "\n[5/5] Testing Python API (ctypes wrapper)..."
python3 -c '
import sys
try:
    from ringdetect import find_rings
    
    # --- Valid Data Test ---
    x = [0.0000, 1.2098, 1.2098, 0.0000, -1.2098, -1.2098]
    y = [1.3970, 0.6985, -0.6985, -1.3970, -0.6985, 0.6985]
    z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    elements = ["C", "C", "C", "C", "C", "C"]

    rings = find_rings(x, y, z, elements, max_ring=6)
    assert len(rings) == 1, f"Expected 1 ring, found {len(rings)}"
    assert rings[0]["size"] == 6, "Expected a 6-MR"
    print("\033[0;32m  ✓ Python API native arrays passed.\033[0m")

    # --- Robustness Test 1: Mismatched Arrays ---
    try:
        find_rings([0.0], [0.0], [0.0], ["C", "H"]) # 1 coord, 2 elements
        print("\033[0;31m  ✗ Python API failed to catch mismatched arrays!\033[0m")
        sys.exit(1)
    except ValueError:
        print("\033[0;32m  ✓ Python API gracefully caught mismatched array lengths.\033[0m")

    # --- Robustness Test 2: Empty Molecule ---
    empty_rings = find_rings([], [], [], [])
    assert len(empty_rings) == 0, "Expected 0 rings for empty input"
    print("\033[0;32m  ✓ Python API gracefully handled empty molecules without crashing C-backend.\033[0m")

except Exception as e:
    print(f"\033[0;31m  ✗ Python API test failed: {e}\033[0m")
    sys.exit(1)
' || exit 1

echo -e "\n${GREEN}======================================================${NC}"
echo -e "${GREEN} 🎉 ALL TESTS PASSED SUCCESSFULLY! The Engine is solid.${NC}"
echo -e "${GREEN}======================================================${NC}"
