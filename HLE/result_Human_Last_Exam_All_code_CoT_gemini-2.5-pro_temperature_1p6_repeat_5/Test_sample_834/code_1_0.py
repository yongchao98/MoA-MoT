import math

# Step 1: Calculate cross-sectional properties
# Given a = 12^(1/4)
# Therefore, a^4 = 12
a_raised_to_4 = 12.0

# The moments of inertia are calculated using the parallel axis theorem.
# I_total = I_main_square - 2 * I_cutout
# I_zz = I about the s-axis, I_ss = I about the z-axis.

# I_zz for the composite section
# Izz = (81*a^4/12) - 2 * (a^4/12 + a^2 * a^2)
# Izz = (81/12)*a^4 - 2 * (13/12)*a^4 = (81 - 26)/12 * a^4 = (55/12) * a^4
I_zz = (55 / 12) * a_raised_to_4

# I_ss for the composite section
# Iss = (81*a^4/12) - 2 * (a^4/12 + a^2 * (a/2)^2)
# Iss = (81/12)*a^4 - 2 * (a^4/12 + a^4/4) = (81/12)*a^4 - 2 * (4/12)*a^4 = (81-8)/12 * a^4 = (73/12) * a^4
I_ss = (73 / 12) * a_raised_to_4

print(f"Calculated moment of inertia I_zz = {I_zz}")
print(f"Calculated moment of inertia I_ss = {I_ss}")
print("-" * 30)

# Step 2: Calculate the given parameters L and q0
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

print(f"Calculated length parameter L = {L}")
print(f"Calculated load parameter q0 = {q0}")
print("-" * 30)

# Step 3 & 4: Derive the expression for F
# The deflection at x = 3L/2 must be zero.
# Total deflection w_total = w_F + w_q = 0
# Upward deflection from F at x=3L/2 is w_F = F * (3L/2)^3 / (3*EI) = 9*F*L^3 / (8*EI)
# Downward deflection from q(x) at x=3L/2 is w_q.
# Deflection at x=L due to q: w_q(L) = -q0*L^4 / (30*EI)
# Slope at x=L due to q: theta_q(L) = -q0*L^3 / (24*EI)
# Extrapolating to x=3L/2: w_q(3L/2) = w_q(L) + theta_q(L) * (L/2)
# w_q(3L/2) = -q0*L^4/(30*EI) - q0*L^4/(48*EI) = - (13/240) * q0*L^4/EI
# Setting w_total = 0:
# 9*F*L^3 / (8*EI) = (13/240) * q0*L^4 / EI
# 9*F / 8 = 13*q0*L / 240
# F = (13 * q0 * L / 240) * (8 / 9)
# F = (13 * q0 * L) / 270

# Step 5: Calculate the final numerical value of F
F_numerator = 13 * q0 * L
F_denominator = 270
F_value = F_numerator / F_denominator

print("The final equation for F is derived as F = (13 * q0 * L) / 270")
print("Substituting the values of q0 and L:")
print(f"F = (13 * {q0} * {L}) / 270")
print(f"F = {F_value}")

print("-" * 30)
print("The final result for F is:")
print(f"<<<{F_value}>>>")