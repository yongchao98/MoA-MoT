import math

# Step 1: Define the crystal parameters for NaMgH3.
# NaMgH3 has a cubic perovskite structure.
# The lattice parameter 'a' at room temperature is approximately 3.69 Å.
a = 3.69

# Step 2: Identify the Miller indices for the second major peak.
# The first major diffraction peak for a perovskite is (110).
# The second major diffraction peak is (200).
h, k, l = 2, 0, 0

# Step 3: Calculate the Q-space position using the standard formula.
# Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)
hkl_sq_sum = h**2 + k**2 + l**2
sqrt_hkl_sq = math.sqrt(hkl_sq_sum)
q_value = (2 * math.pi / a) * sqrt_hkl_sq

# Step 4: Print the explanation and the final equation with all numbers.
print(f"The calculation for the Q-space position of the ({h},{k},{l}) peak is:")
print(f"Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)")
print("\nSubstituting the known values:")
print(f"Lattice parameter 'a' = {a} Å")
print(f"Miller indices (h, k, l) = ({h}, {k}, {l})\n")
print("The final equation is:")
print(f"Q = (2 * {math.pi:.5f} / {a}) * sqrt({h}^2 + {k}^2 + {l}^2)")
print(f"Q = ({2 * math.pi / a:.5f}) * sqrt({hkl_sq_sum})")
print(f"Q = ({2 * math.pi / a:.5f}) * {sqrt_hkl_sq:.5f}")
print(f"Q = {q_value:.4f} 1/Å")