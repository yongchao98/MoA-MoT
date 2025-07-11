import math

# --- Plan ---
# 1. Define the crystal parameters for NaMgH3 (orthorhombic, Pnma).
# 2. Identify the Miller indices (h,k,l) for the second major diffraction peak.
#    Based on known diffraction data, this corresponds to the doublet of (121) and (002).
# 3. Calculate the Q-space position for one of these reflections, (0,0,2), which represents the location of the peak feature.

# --- Crystal parameters for NaMgH3 ---
a = 5.760  # Angstrom
b = 8.058  # Angstrom
c = 5.670  # Angstrom

# Miller indices for one of the reflections in the second major peak doublet
h, k, l = (0, 0, 2)

print(f"To find the Q-space position of the second major diffraction peak of NaMgH3, we first identify the corresponding Miller indices (h,k,l).")
print("NaMgH3 has an orthorhombic structure (Pnma). Based on its calculated diffraction pattern, the second most intense peak is a doublet composed of the (1,2,1) and (0,0,2) reflections.")
print("We will calculate the Q-space position for the (0,0,2) reflection as a representative for this peak.\n")

# --- Calculation ---

# 1. Calculate the inverse square of the interplanar spacing (d)
inv_d_squared = (h / a)**2 + (k / b)**2 + (l / c)**2

# 2. Calculate the d-spacing
d_spacing = 1 / math.sqrt(inv_d_squared)

# 3. Calculate the scattering vector magnitude (Q)
q_value = (2 * math.pi) / d_spacing

# --- Output Results ---

print("The calculation steps are as follows:")
print(f"1. Calculate 1/d^2 for (h,k,l) = ({h},{k},{l}):")
print(f"   1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2")
print(f"   1/d^2 = ({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2")
print(f"   1/d^2 = {inv_d_squared:.6f} Å⁻²\n")

print(f"2. Calculate the d-spacing:")
print(f"   d = 1 / sqrt(1/d^2)")
print(f"   d = 1 / sqrt({inv_d_squared:.6f}) = {d_spacing:.4f} Å\n")

print(f"3. Calculate the Q-space position:")
print(f"   Q = 2 * pi / d")
print(f"   Q = 2 * {math.pi:.4f} / {d_spacing:.4f} = {q_value:.4f} Å⁻¹")