import math

# Step 1: Define physical constants and standard properties for a monoclonal antibody (mAb).
# These are standard assumptions for a typical IgG-type mAb.
Na = 6.022e23  # Avogadro's number, in mol^-1
M = 150000.0   # Molar mass of a typical mAb, in g/mol
v_bar = 0.73   # Partial specific volume of a typical protein, in mL/g (same as cm^3/g)

# The total experimental B22 of -7.585 mL/g is context, showing that attractive forces
# dominate in this specific condition. The steric component, however, is always positive
# and is calculated based on the molecule's physical dimensions.

print("--- Calculation of the Steric-Only Second Osmotic Virial Coefficient (B22_steric) ---")
print("\nStep 1: Define constants and assumed mAb properties.")
print(f"Avogadro's number (Na): {Na:.3e} mol^-1")
print(f"Molar Mass (M): {M:.1f} g/mol")
print(f"Partial Specific Volume (v_bar): {v_bar:.2f} mL/g")

# Step 2: Calculate the effective radius (r) of the mAb modeled as a hard sphere.
# The formula is derived from M * v_bar = Na * (4/3)*pi*r^3.
r_cubed_cm3 = (3 * M * v_bar) / (4 * math.pi * Na)
r_cm = r_cubed_cm3**(1/3)

print("\nStep 2: Calculate the effective radius (r) of the mAb.")
print("Equation for r^3: (3 * M * v_bar) / (4 * pi * Na)")
print(f"r^3 = (3 * {M:.1f} * {v_bar:.2f}) / (4 * {math.pi:.4f} * {Na:.3e})")
print(f"r^3 = {r_cubed_cm3:.4e} cm^3")
print(f"r = {r_cm:.4e} cm")

# Step 3: Calculate B22_steric using the hard-sphere model.
# The formula gives the result in mL/g, as required.
b22_steric_numerator = 16 * math.pi * r_cubed_cm3 * Na
b22_steric_denominator = 3 * M
b22_steric = b22_steric_numerator / b22_steric_denominator

print("\nStep 3: Calculate B22_steric.")
print("Equation: (16 * pi * r^3 * Na) / (3 * M)")
print(f"B22_steric = (16 * {math.pi:.4f} * {r_cubed_cm3:.4e} * {Na:.3e}) / (3 * {M:.1f})")
print(f"B22_steric = {b22_steric_numerator:.4e} / {b22_steric_denominator:.1f}")

print("\n--- Final Answer ---")
print(f"The second osmotic virial coefficient from Steric-only behavior for the monoclonal antibody is {b22_steric:.3f} mL/g.")

# For clarity, let's break down the final equation with all numbers:
# B22_steric = (16 * 3.1416 * 4.3413e-20 cm^3 * 6.022e+23 mol^-1) / (3 * 150000.0 g/mol)
# B22_steric = 1313437.5 cm^3/mol / 450000.0 g/mol
# B22_steric = 2.919 cm^3/g, which is 2.919 mL/g.
