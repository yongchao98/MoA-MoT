import numpy as np

# Step 1: Define constants and calculate geometric properties
# a = 12^(1/4), so a^4 = 12
a_4 = 12.0
a = np.power(a_4, 1/4)

# Moment of inertia for the large square (3a x 3a) about its centroid (origin)
I_ss_large = (3 * a) * (3 * a)**3 / 12
I_zz_large = (3 * a) * (3 * a)**3 / 12

# Moment of inertia for one small square cutout (a x a) using Parallel Axis Theorem
# I = I_centroid + A * d^2
area_small = a**2
I_ss_small_centroid = a**4 / 12
I_zz_small_centroid = a**4 / 12

# For the cutouts at (s=a/2, z=-a) and (s=-a/2, z=a), the distances from the main axes are:
d_s = a / 2
d_z = a
# By symmetry, the contribution from both cutouts is identical.
I_ss_cutouts = 2 * (I_ss_small_centroid + area_small * d_z**2)
I_zz_cutouts = 2 * (I_zz_small_centroid + area_small * d_s**2)

# Total moments of inertia for the composite cross-section
I_ss = I_ss_large - I_ss_cutouts
I_zz = I_zz_large - I_zz_cutouts

# Step 2: Calculate L and q0
# Using the calculated I_ss and I_zz to find L and q0
L = 30 * I_zz / 73
q0 = 9 * I_ss / 55

# Step 3: Use beam deflection formulas to find the relation for F
# Deflection at tip due to triangular load q:
# y_q(L) = -q0 * L**4 / (30 * EI)
# theta_q(L) = -q0 * L**3 / (24 * EI)
# y_q(3L/2) = y_q(L) + theta_q(L) * (L/2) = (-q0*L**4/EI)*(1/30 + 1/48) = -13*q0*L**4/(240*EI)
# Deflection at tip due to point force F:
# y_F(3L/2) = F * (3*L/2)**3 / (3*EI) = 9*F*L**3/(8*EI)
# Setting y_q + y_F = 0 gives:
# 9*F*L**3 / 8 = 13*q0*L**4 / 240
# F = (13 * q0 * L / 240) * (8 / 9) = 13 * q0 * L / 270

# Step 4: Calculate the final numerical value for F
F = (13 * q0 * L) / 270

# Step 5: Print the results including the final equation
print(f"Based on the geometry, the moments of inertia are I_ss = {I_ss:.2f} and I_zz = {I_zz:.2f}.")
print(f"This gives the parameters L = {L:.2f} and q0 = {q0:.2f}.")
print("-" * 30)
print("The final equation for the force F is derived as:")
print("F = 13 * q0 * L / 270")
print("-" * 30)
print("Substituting the values:")
print(f"F = 13 * {q0:.2f} * {L:.2f} / 270")
print(f"F = {F:.2f}")

# Final Answer
print(f"\n<<<{F:.2f}>>>")