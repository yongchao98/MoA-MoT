import math

# Step 1 & 2: Derive the equation for F
# The upward deflection at the tip (x=3L/2) due to F is y_F = F * (3L/2)^3 / (3*EI) = 9*F*L^3 / (8*EI).
# The downward deflection at the tip due to the triangular load q(x) is y_q(3L/2) = y_q(L) + theta_q(L)*(L/2).
# Using standard formulas: y_q(L) = q_0*L^4/(30*EI) and theta_q(L) = q_0*L^3/(24*EI).
# y_q(3L/2) = q_0*L^4/(30*EI) + q_0*L^3/(24*EI)*(L/2) = q_0*L^4/(30*EI) + q_0*L^4/(48*EI)
# y_q(3L/2) = (1/30 + 1/48) * q_0*L^4/EI = ((16+10)/480) * q_0*L^4/EI = 26/480 * q_0*L^4/EI = 13*q_0*L^4/(240*EI).
# Set deflections equal: 9*F*L^3/(8*EI) = 13*q_0*L^4/(240*EI).
# Solving for F: F = (13 * q_0 * L / 240) * (8 / 9) = 13 * q_0 * L / 270.

# Step 3: Calculate geometric properties I_ss and I_zz
# Given a = 12^(1/4), so a^4 = 12
a_fourth = 12.0

# Calculate I_zz = I_zz(large_square) - 2 * I_zz(cutout)
# I_zz_large = (3a)^4 / 12 = 81*a^4 / 12
# I_zz_cutout = I_c + A*d^2 = a^4/12 + a^2*a^2 = 13*a^4/12
# I_zz = (81*a^4/12) - 2*(13*a^4/12) = (81-26)*a^4/12 = 55*a^4/12
I_zz = 55 * a_fourth / 12
print(f"Calculated I_zz = {I_zz}")

# Calculate I_ss = I_ss(large_square) - 2 * I_ss(cutout)
# I_ss_large = (3a)^4 / 12 = 81*a^4 / 12
# I_ss_cutout = I_c + A*d^2 = a^4/12 + a^2*(a/2)^2 = a^4/12 + a^4/4 = 4*a^4/12
# I_ss = (81*a^4/12) - 2*(4*a^4/12) = (81-8)*a^4/12 = 73*a^4/12
I_ss = 73 * a_fourth / 12
print(f"Calculated I_ss = {I_ss}")
print("-" * 20)

# Step 4: Calculate L and q_0
# Given L = 30*I_zz/73 and q_0 = 9*I_ss/55
L = 30 * I_zz / 73
q_0 = 9 * I_ss / 55

print(f"Calculated L = {L}")
print(f"Calculated q_0 = {q_0}")
print("-" * 20)

# Step 5: Solve for F
# F = 13 * q_0 * L / 270
F = (13 / 270) * q_0 * L

print("The final equation for F is:")
print(f"F = (13 / 270) * q_0 * L")
print(f"F = (13 / 270) * ({q_0:.4f}) * ({L:.4f})")
print(f"F = {F}")
print("-" * 20)
print(f"The required force F is {F}.")
print("<<<13.0>>>")