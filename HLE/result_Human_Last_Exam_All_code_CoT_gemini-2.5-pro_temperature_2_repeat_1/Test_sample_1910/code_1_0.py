import math

# Step 1: Define crystal structure parameters for orthorhombic NaMgH3 at room temperature.
a = 5.93  # Lattice parameter a in Angstroms
b = 8.32  # Lattice parameter b in Angstroms
c = 5.88  # Lattice parameter c in Angstroms

# Step 2: Define Miller indices for a representative reflection of the second major diffraction peak.
# The second major peak is a multiplet composed of the (220) and (022) reflections.
# We will calculate the position for the (022) plane.
h, k, l = 0, 2, 2

print("Calculation for the Q-space position of the second major diffraction peak of NaMgH3.")
print(f"The crystal structure is orthorhombic with lattice parameters: a = {a} Å, b = {b} Å, c = {c} Å.")
print(f"We are calculating the position for the (hkl) = ({h}{k}{l}) reflection, which is part of the second major peak group.")
print("-" * 30)

# Step 3: Calculate the interplanar spacing 'd' for the (022) plane.
print("First, we calculate the interplanar spacing 'd'.")
print("The formula for an orthorhombic system is: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2")
# Calculation
h_a_sq = (h/a)**2
k_b_sq = (k/b)**2
l_c_sq = (l/c)**2
d_squared_inv = h_a_sq + k_b_sq + l_c_sq
d = math.sqrt(1 / d_squared_inv)
# Outputting the numbers for the formula
print(f"1/d^2 = ({h}/{a})^2 + ({k}/{b})^2 + ({l}/{c})^2")
print(f"1/d^2 = {h_a_sq:.4f} + {k_b_sq:.4f} + {l_c_sq:.4f}")
print(f"1/d^2 = {d_squared_inv:.4f} Å^-2")
print(f"d = 1 / sqrt({d_squared_inv:.4f}) = {d:.4f} Å")
print("-" * 30)

# Step 4: Calculate the Q-space position.
print("Next, we calculate the Q-space position using the formula: Q = 2 * pi / d")
pi = math.pi
Q = (2 * pi) / d
# Outputting the numbers for the formula
print(f"Q = 2 * {pi:.4f} / {d:.4f}")
print(f"The Q-space position for the ({h}{k}{l}) reflection is: {Q:.4f} Å^-1.")
