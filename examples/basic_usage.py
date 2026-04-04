from ringdetect.engine import find_rings

# Pass arrays directly (no file reading necessary!)
x = [0.0000, 1.2098, 1.2098, 0.0000, -1.2098, -1.2098]
y = [1.3970, 0.6985, -0.6985, -1.3970, -0.6985, 0.6985]
z = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
elements = ["C", "C", "C", "C", "C", "C"]

rings = find_rings(x, y, z, elements, max_ring=6, threads=4)

for ring in rings:
    print(ring)
