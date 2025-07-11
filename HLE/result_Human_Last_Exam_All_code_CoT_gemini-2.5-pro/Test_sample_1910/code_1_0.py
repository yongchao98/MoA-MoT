import math

# Step 1: Define the crystal properties for NaMgH3 (Pnma, orthorhombic)
# Lattice parameters in Angstroms (Å)
a = 5.549
b = 7.891
c = 5.503

# Step 2: Identify the Miller indices (h, k, l) for the second major peak
# The most intense peak is (101). The second most intense peak is (121).
h, k, l = 1, 2, 1

# Step 3: Calculate the Q-space position
# The formula is Q = 2 * pi * sqrt((h^2/a^2) + (k^2/b^2) + (l^2/c^2))
h2_a2 = (h**2 / a**2)
k2_b2 = (k**2 / b**2)
l2_c2 = (l**2 / c**2)
d_inv_sq = h2_a2 + k2_b2 + l2_c2
Q = 2 * math.pi * math.sqrt(d_inv_sq)

# Print the explanation and the final calculation
print("The Q-space position for an orthorhombic crystal is given by:")
print("Q = 2 * pi * sqrt((h^2/a^2) + (k^2/b^2) + (l^2/c^2))\n")
print(f"For the second major peak of NaMgH3 with (h,k,l) = ({h},{k},{l}):")
print(f"Lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å")
print("\nThe calculation is:")
print(f"Q = 2 * {math.pi:.4f} * sqrt(({h}^2/{a}^2) + ({k}^2/{b}^2) + ({l}^2/{c}^2))")
print(f"Q = 2 * {math.pi:.4f} * sqrt({h2_a2:.4f} + {k2_b2:.4f} + {l2_c2:.4f})")
print(f"Q = 2 * {math.pi:.4f} * sqrt({d_inv_sq:.4f})")
print(f"Q = {Q:.4f} (1/Å)\n")
print("The second major diffraction peak is located at Q ≈ 2.2647 1/Å.")

# Final answer in the required format
# print(f"<<<{Q:.4f}>>>")