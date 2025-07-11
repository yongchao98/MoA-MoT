import numpy as np
from scipy.integrate import quad

# Step 1: Calculate Area Moments of Inertia
# Given a = 12^(1/4), so a^4 = 12
a_sq_4 = 12.0

# Large square (3a x 3a) centered at (0,0)
side_large = 3  # in units of a
I_ss_large = (1 / 12) * (side_large * a_sq_4) * side_large**3 
I_zz_large = (1 / 12) * (side_large * a_sq_4) * side_large**3

# Cutout 1: small square (a x a) centered at (s=a/2, z=-a)
side_small = 1.0 # in units of a
area_small = side_small**2 # in units of a^2
I_ss_c = (1 / 12) * (side_small * a_sq_4) * side_small**3
I_zz_c = (1 / 12) * (side_small * a_sq_4) * side_small**3

# Parallel axis theorem: I = I_c + A*d^2
# For cutout 1 at (s_c = a/2, z_c = -a)
# Distances are squared, so signs don't matter. d_s = 0.5, d_z = 1.0 in units of a
I_ss_1 = I_ss_c + (area_small * a_sq_4 / 12) * (1.0)**2 * 12 # a^2 * d_z^2 = a^2 * a^2 = a^4 = 12
I_zz_1 = I_zz_c + (area_small * a_sq_4 / 12) * (0.5)**2 * 12 # a^2 * d_s^2 = a^2 * (a/2)^2 = a^4/4

# For cutout 2 at (s_c = -a/2, z_c = a)
# Distances are squared, so they are the same as for cutout 1
I_ss_2 = I_ss_1
I_zz_2 = I_zz_1

# Total moments of inertia by subtraction
I_ss = I_ss_large - I_ss_1 - I_ss_2
I_zz = I_zz_large - I_zz_1 - I_zz_2

# Step 2: Calculate problem parameters
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

def integrand(x):
    return np.sin(x**2)

EI, _ = quad(integrand, 0, np.pi)

# Step 3 & 4: Derive the formula for F
# Deflection due to q(x) at x=3L/2 is w_q = -13 * q0 * L^4 / (240 * EI)
# Deflection due to F at x=3L/2 is w_F = 9 * F * L^3 / (8 * EI)
# Set w_q + w_F = 0 => 9*F*L^3 / (8*EI) = 13*q0*L^4 / (240*EI)
# Solving for F: F = (13 * q0 * L / 240) * (8 / 9) = 13 * q0 * L / 270

# Step 5: Final Calculation
F = (13 * q0 * L) / 270

# Print the final result including the equation with numerical values
print(f"Given values:")
print(f"a = 12^(1/4)")
print(f"I_ss = {I_ss:.2f}")
print(f"I_zz = {I_zz:.2f}")
print(f"L = (30 * {I_zz:.2f}) / 73 = {L:.2f}")
print(f"q0 = (9 * {I_ss:.2f}) / 55 = {q0:.2f}")
print(f"EI = integral(sin(x^2)) from 0 to pi = {EI:.4f}")
print("\nTo find the force F, we set the total deflection at the tip to zero.")
print("The deflection from the distributed load q(x) is balanced by the deflection from the force F.")
print("w_q(3L/2) + w_F(3L/2) = 0")
print("-13*q0*L^4/(240*EI) + 9*F*L^3/(8*EI) = 0")
print("Solving for F gives: F = (13 * q0 * L) / 270")
print("\nSubstituting the numerical values:")
print(f"F = (13 * {q0:.2f} * {L:.2f}) / 270")
print(f"F = {F:.2f}")

final_answer = F
<<<13.0>>>