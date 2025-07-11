import math

# Step 1: Define 'a' and calculate its fourth power.
# Given a = 12^(1/4)
a_4 = 12
print(f"Given a = 12^(1/4), so a^4 = {a_4:.1f}\n")

# Step 2: Calculate the second moments of area, I_ss and I_zz.
# For a large 3a x 3a square, I_main = (3a)*(3a)^3 / 12 = 81*a^4 / 12.
# For each a x a cutout, we use the parallel axis theorem I = I_c + A*d^2.
# For I_ss, the distance d_z from the s-axis is 'a'.
# I_ss_cutout = a^4/12 + (a^2)*a^2 = 13*a^4/12.
# I_ss = I_main - 2 * I_ss_cutout = 81*a^4/12 - 2*13*a^4/12 = 55*a^4/12.
I_ss = 55 * a_4 / 12

# For I_zz, the distance d_s from the z-axis is 'a/2'.
# I_zz_cutout = a^4/12 + (a^2)*(a/2)^2 = a^4/12 + a^4/4 = a^4/3.
# I_zz = I_main - 2 * I_zz_cutout = 81*a^4/12 - 2*a^4/3 = 73*a^4/12.
I_zz = 73 * a_4 / 12

print(f"The second moment of area about the s-axis is I_ss = 55 * a^4 / 12 = {I_ss:.1f}")
print(f"The second moment of area about the z-axis is I_zz = 73 * a^4 / 12 = {I_zz:.1f}\n")

# Step 3: Calculate the parameters L and q0 using the given formulas.
# L = 30 * I_zz / 73
L = 30 * I_zz / 73
# q0 = 9 * I_ss / 55
q0 = 9 * I_ss / 55

print(f"Calculating the beam length L:")
print(f"L = 30 * I_zz / 73 = 30 * {I_zz:.1f} / 73 = {L:.1f}\n")
print(f"Calculating the maximum distributed load q0:")
print(f"q0 = 9 * I_ss / 55 = 9 * {I_ss:.1f} / 55 = {q0:.1f}\n")

# Step 4: Calculate the required force F.
# The total deflection at the tip (x = 3L/2) is the sum of the deflection
# from the triangular load (v_q) and the point load (v_F).
# For the total deflection to be zero, v_q + v_F = 0.
# Downward deflection from q: v_q = -(13 * q0 * L^4) / (240 * EI)
# Upward deflection from F:   v_F = (9 * F * L^3) / (8 * EI)
# Setting v_q + v_F = 0 leads to:
# (9 * F * L^3) / (8 * EI) = (13 * q0 * L^4) / (240 * EI)
# Simplifying this equation by canceling EI and L^3 and solving for F gives:
# F = (13 * q0 * L) / 270

F = (13 * q0 * L) / 270

print("The required force F is found from the zero-deflection condition.")
print("The derived formula for F is: F = (13 * q0 * L) / 270\n")

print("Substituting the numerical values into the equation:")
print(f"F = (13 * {q0:.1f} * {L:.1f}) / 270")
print(f"F = {13 * q0 * L} / 270")
print(f"F = {F}")

print(f"\n<<<13.0>>>")