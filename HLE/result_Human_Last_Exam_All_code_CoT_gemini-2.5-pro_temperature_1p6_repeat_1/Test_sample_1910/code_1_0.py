import math

# Step 1: Define the parameters for NaMgH3 perovskite using a cubic approximation.
# The lattice parameter 'a' for the pseudo-cubic cell of NaMgH3 is approximately 3.89 Å.
a = 3.89  # units: Angstroms

# Step 2: Identify the Miller indices (h, k, l) for the second major diffraction peak.
# For a simple cubic perovskite structure, the peaks appear in order of increasing h²+k²+l².
# 1st peak: (100) -> h²+k²+l² = 1
# 2nd peak: (110) -> h²+k²+l² = 2
# Therefore, we use the Miller indices for the (110) peak.
h, k, l = 1, 1, 0

# Step 3: Calculate the sum of the squares of the Miller indices.
hkl_sq_sum = h**2 + k**2 + l**2

# Step 4: Calculate the d-spacing for the (110) plane using the formula for a cubic lattice.
# The formula is: d = a / sqrt(h² + k² + l²)
try:
    d_spacing = a / math.sqrt(hkl_sq_sum)
except ZeroDivisionError:
    d_spacing = float('inf')


# Step 5: Calculate the Q-space position.
# The magnitude of the scattering vector Q is related to the d-spacing by Q = 2π / d.
# The provided wavelength is not needed for this direct calculation.
if d_spacing != float('inf'):
    Q = 2 * math.pi / d_spacing
else:
    Q = 0

# Step 6: Print the steps and the final result.
print("--- Calculation for the Second Diffraction Peak of NaMgH3 ---")
print(f"Approximated Crystal Structure: Cubic Perovskite")
print(f"Lattice parameter (a): {a} Å")
print(f"Miller indices (h,k,l) for the second peak: ({h}, {k}, {l})")
print("\n--- Intermediate Calculations ---")
print(f"Equation for d-spacing: d = a / sqrt(h² + k² + l²)")
print(f"d = {a} / sqrt({h}² + {k}² + {l}²) = {d_spacing:.4f} Å")
print("\n--- Final Q-space Calculation ---")
print(f"Equation for Q: Q = 2 * pi / d")
print(f"Q = 2 * {math.pi:.4f} / {d_spacing:.4f}")
print(f"\nThe position of the second major diffraction peak in Q-space is: {Q:.4f} 1/Å")