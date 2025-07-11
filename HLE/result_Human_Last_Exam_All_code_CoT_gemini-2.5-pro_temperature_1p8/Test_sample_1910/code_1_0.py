import math

# Step 1: Define the properties of NaMgH3 perovskite.
# Lattice parameter 'a' for cubic NaMgH3 at room temperature is ~3.65 Angstroms.
lattice_parameter = 3.65  # units: Angstroms

# Step 2: Identify the Miller indices (h, k, l) for the second major peak.
# The first major peak is (110). The second major peak is (200).
h, k, l = 2, 0, 0

# Step 3: Calculate the Q-space position.
# The formula is Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2).
hkl_sum_sq = h**2 + k**2 + l**2
q_value = (2 * math.pi / lattice_parameter) * math.sqrt(hkl_sum_sq)

# Step 4: Print the final calculation and the result.
# The problem asks to output each number in the final equation.
print(f"The second major diffraction peak for NaMgH3 corresponds to the Miller indices (h, k, l) = ({h}, {k}, {l}).")
print("The Q-space position is calculated using the formula: Q = (2 * pi / a) * sqrt(h^2 + k^2 + l^2)")
print("\nPlugging in the values:")
print(f"Q = (2 * {math.pi:.5f} / {lattice_parameter}) * sqrt({h}^2 + {k}^2 + {l}^2)")
print(f"Q = ({2 * math.pi:.5f} / {lattice_parameter}) * sqrt({hkl_sum_sq})")
print(f"Q = {2 * math.pi / lattice_parameter:.5f} * {math.sqrt(hkl_sum_sq):.5f}")
print(f"\nThe calculated Q-space position is: {q_value:.4f} 1/Angstrom")