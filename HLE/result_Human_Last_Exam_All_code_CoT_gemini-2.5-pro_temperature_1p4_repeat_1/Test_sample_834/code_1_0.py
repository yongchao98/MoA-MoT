import numpy as np

# Step 1: Calculate the moments of inertia
# Given a = 12^(1/4), so a^4 = 12
a_sq_4 = 12.0

# I_zz = integral(s^2 dA)
# For a rectangle with base h (s-dim) and height b (z-dim), Izz_c = b * h^3 / 12
# For a rectangle with base h (s-dim) and height b (z-dim), I_ss_c = h * b^3 / 12

# Main square (3a x 3a) centered at origin
I_zz_large_sq = (3 * np.sqrt(a_sq_4)) * (3 * np.sqrt(a_sq_4))**3 / 12
I_ss_large_sq = (3 * np.sqrt(a_sq_4))**3 * (3 * np.sqrt(a_sq_4)) / 12
# simplified:
I_zz_large_sq = 81 * a_sq_4 / 12
I_ss_large_sq = 81 * a_sq_4 / 12

# Cutout square 1 (a x a) centered at (s=a/2, z=-a)
# Centroidal moment of inertia
I_zz_c_small = a_sq_4 / 12
I_ss_c_small = a_sq_4 / 12
# Area
A_small = np.sqrt(a_sq_4) ** 2
# Parallel axis theorem: I = I_c + A * d^2
# For cutout 1, d_s = a/2, d_z = -a
a = np.sqrt(a_sq_4)
I_zz_cut1 = I_zz_c_small + A_small * (a / 2)**2
I_ss_cut1 = I_ss_c_small + A_small * (-a)**2

# Cutout square 2 (a x a) centered at (s=-a/2, z=a)
# For cutout 2, d_s = -a/2, d_z = a
I_zz_cut2 = I_zz_c_small + A_small * (-a / 2)**2
I_ss_cut2 = I_ss_c_small + A_small * (a)**2

# Total moments of inertia
I_zz = I_zz_large_sq - I_zz_cut1 - I_zz_cut2
I_ss = I_ss_large_sq - I_ss_cut1 - I_ss_cut2

print(f"Moments of Inertia Calculation:")
print(f"a^4 = {a_sq_4}")
print(f"I_zz for the composite section = {I_zz:.2f}")
print(f"I_ss for the composite section = {I_ss:.2f}\n")


# Step 2: Calculate L and q_0
L = 30 * I_zz / 73
q_0 = 9 * I_ss / 55

print("Given Data Calculation:")
print(f"L = (30 * {I_zz:.2f}) / 73 = {L:.2f}")
print(f"q_0 = (9 * {I_ss:.2f}) / 55 = {q_0:.2f}\n")


# Step 3 & 4: Derive and solve for F
# The condition for zero deflection at x=3L/2 is w_F(3L/2) + w_q(3L/2) = 0
# w_F(3L/2) is the deflection due to force F, which is upwards (positive)
# w_q(3L/2) is the deflection due to the distributed load q, which is downwards (negative)
# w_F(3L/2) = (9 * F * L^3) / (8 * EI)
# w_q(3L/2) = (-37 * q_0 * L^4) / (240 * EI)
# Setting sum to zero: (9 * F * L^3) / (8 * EI) = (37 * q_0 * L^4) / (240 * EI)
# Solving for F: F = (37 * q_0 * L) / 270

F = (37 * q_0 * L) / 270

print("Force Calculation:")
print("The relation for F is derived from setting the total deflection at the tip to zero:")
print("Deflection(F) + Deflection(q) = 0")
print("(9 * F * L^3) / (8 * EI) - (37 * q_0 * L^4) / (240 * EI) = 0")
print("Solving for F gives: F = (37 * q_0 * L) / 270\n")
print("Substituting the values:")
print(f"F = (37 * {q_0:.2f} * {L:.2f}) / 270")
print(f"F = {37 * q_0 * L:.2f} / 270")
print(f"F = {F:.2f}")

final_answer = F
print(f"\nThe required force F for zero deflection at the endpoint is {final_answer:.2f}.")
print(f"<<<{final_answer}>>>")