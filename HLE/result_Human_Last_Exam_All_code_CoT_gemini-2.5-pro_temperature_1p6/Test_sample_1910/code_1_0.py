import math

# Step 1: Define the crystal structure parameters for orthorhombic NaMgH3
a = 5.530  # lattice parameter a in Angstrom
b = 7.811  # lattice parameter b in Angstrom
c = 5.518  # lattice parameter c in Angstrom

# Step 2: Define Miller indices for the second major diffraction peak
# Based on diffraction patterns for this structure, the second major peak is (202).
h, k, l = 2, 0, 2

print(f"Calculating the Q-space position for the perovskite NaMgH3.")
print(f"Structure: Orthorhombic, with lattice parameters a={a} Å, b={b} Å, c={c} Å.")
print(f"The second major diffraction peak corresponds to the Miller indices (h,k,l) = ({h},{k},{l}).")
print("The Q-space position is independent of the measurement wavelength.\n")

# Step 3: Calculate the interplanar spacing (d)
# Formula for orthorhombic system: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2
d_sq_inv = (h / a)**2 + (k / b)**2 + (l / c)**2
d_spacing = math.sqrt(1 / d_sq_inv)

print(f"Equation for d-spacing: 1/d^2 = ({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2")
print(f"The calculated interplanar d-spacing for the ({h},{k},{l}) plane is: {d_spacing:.4f} Å")

# Step 4: Calculate the Q-space position
# Formula: Q = 2 * pi / d
q_value = (2 * math.pi) / d_spacing

print(f"\nEquation for Q-space position: Q = 2 * pi / d")
print(f"Final calculation: Q = (2 * {math.pi:.4f}) / {d_spacing:.4f}")
print(f"\nThe Q-space position for the second major diffraction peak is: {q_value:.4f} 1/Å")