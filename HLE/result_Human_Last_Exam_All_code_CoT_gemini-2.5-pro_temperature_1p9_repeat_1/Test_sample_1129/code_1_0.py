import math

# Step 1: Define the constants and given values based on the problem description and plan.

# The Lawson criterion for ignition (n*tau), using the optimistic value for D-T fusion.
# Unit: s * m^-3
lawson_product = 1e20

# Given particle confinement time (tau).
# Unit: s
confinement_time = 1.0

# Given side length of the cubic reaction chamber.
# Unit: cm
chamber_side_length_cm = 10.0

# Step 2: Convert units to be consistent (meters for length).
# Convert side length from cm to m.
chamber_side_length_m = chamber_side_length_cm / 100.0

# Step 3: Calculate the required particle density (n) from the Lawson criterion.
# n = (n * tau) / tau
required_density = lawson_product / confinement_time

# Step 4: Calculate the volume of the reaction chamber.
# V = L^3
chamber_volume = chamber_side_length_m ** 3

# Step 5: Calculate the total minimum number of atoms required.
# N = n * V
total_atoms = required_density * chamber_volume

# Step 6: Print the final result in a clear, descriptive format,
# showing the equation and the values used.
print("This calculation determines the minimum number of atoms for a hypothetical sustained fusion reaction.")
print("The calculation is based on the Lawson criterion, under the optimistic conditions specified.")
print("\n--- Calculation Steps ---")
print(f"1. Minimum Lawson Criterion Product (n\u03C4): {lawson_product:.1e} s/m\u00b3")
print(f"2. Given Confinement Time (\u03C4): {confinement_time} s")
print(f"3. Given Chamber Side Length (L): {chamber_side_length_cm} cm = {chamber_side_length_m} m")

print("\n--- Final Equation ---")
# Using scientific notation for clarity.
print(f"Minimum Atoms (N) = (Lawson Product / Confinement Time) \u00d7 (Side Length)\u00b3")
print(f"N = ({lawson_product:.1e} s/m\u00b3 / {confinement_time} s) \u00d7 ({chamber_side_length_m} m)\u00b3")
print(f"N = {required_density:.1e} m\u207b\u00b3 \u00d7 {chamber_volume:.1e} m\u00b3")
print(f"\nMinimum number of Ti-50 atoms required = {total_atoms:.1e}")
