import numpy as np

# Step 1: Define the given value for 'a'
a = 12**(1/4)
a_p4 = 12  # a^4 = 12

# Step 2: Calculate the second moments of area I_ss and I_zz
# The cross-section is a large 3a x 3a square with two a x a squares cut out.
# We use the parallel axis theorem: I = I_centroid + A*d^2

# Properties of the large square (3a x 3a) centered at the origin
I_ss_large = (3*a) * (3*a)**3 / 12 = 27 * a_p4 / 4
I_zz_large = (3*a)**3 * (3*a) / 12 = 27 * a_p4 / 4

# Properties of the two cut-outs (a x a)
# Cut-out 1 centered at (s=a/2, z=-a)
A_cutout = a**2
I_c_cutout = a**4 / 12
# For I_ss, distance d_z = -a
I_ss_c1 = I_c_cutout + A_cutout * (-a)**2 = a_p4 / 12 + a**2 * a**2 = 13 * a_p4 / 12
# For I_zz, distance d_s = a/2
I_zz_c1 = I_c_cutout + A_cutout * (a/2)**2 = a_p4 / 12 + a**2 * (a**2 / 4) = a_p4 / 3

# Cut-out 2 centered at (s=-a/2, z=a)
# For I_ss, distance d_z = a
I_ss_c2 = I_c_cutout + A_cutout * (a)**2 = 13 * a_p4 / 12
# For I_zz, distance d_s = -a/2
I_zz_c2 = I_c_cutout + A_cutout * (-a/2)**2 = a_p4 / 3

# Total moments of inertia by subtracting the cut-outs
I_ss = I_ss_large - I_ss_c1 - I_ss_c2
I_zz = I_zz_large - I_zz_c1 - I_zz_c2

# Step 3: Calculate the numerical values for L and q0
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

# Step 4: Use the derived formula to find F
# The condition for zero deflection at x = 3L/2 is:
# v_q(3L/2) + v_F(3L/2) = 0
# -(13 * q0 * L^4) / (240 * EI) + (9 * F * L^3) / (8 * EI) = 0
# Solving for F gives: F = (13 * q0 * L) / 270

c1 = 13
c2 = 270
F = (c1 * q0 * L) / c2

# Step 5: Print the results as requested
print("The final equation for the force F is: F = (c1 * q0 * L) / c2")
print("The numbers in the final equation are:")
print(f"c1 = {c1}")
print(f"q0 = {q0}")
print(f"L = {L}")
print(f"c2 = {c2}")
print("\nThe calculated force F is:")
print(F)

# Final answer in the specified format
print(f"\n<<<{F}>>>")