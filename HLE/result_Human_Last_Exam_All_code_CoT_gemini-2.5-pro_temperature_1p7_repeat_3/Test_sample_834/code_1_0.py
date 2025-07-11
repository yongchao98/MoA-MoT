import math

# Step 1: Calculate the value of a^4
# a = 12^(1/4)
a_4 = 12

# Step 2: Calculate the moments of inertia I_zz and I_ss
# Using the parallel axis theorem for the hollow cross-section.
# I_zz_total = I_zz_large_square - 2 * I_zz_small_square
# I_zz = (81/12)*a^4 - 2 * (a^4/12 + a^2 * (a/2)^2) = (73/12)*a^4
I_zz = (73/12) * a_4

# I_ss_total = I_ss_large_square - 2 * I_ss_small_square
# I_ss = (81/12)*a^4 - 2 * (a^4/12 + a^2 * a^2) = (55/12)*a^4
I_ss = (55/12) * a_4

# Step 3: Calculate the beam length L and the maximum load q_0
L = (30 * I_zz) / 73
q_0 = (9 * I_ss) / 55

# Step 4: Calculate the required force F using the derived formula
# F = (13 * q_0 * L) / 270
F = (13 * q_0 * L) / 270

# Step 5: Print the results, including the final equation with numerical values
print(f"Calculated value for a^4: {a_4}")
print(f"Moment of inertia I_zz: {I_zz}")
print(f"Moment of inertia I_ss: {I_ss}")
print(f"Beam length L: {L}")
print(f"Maximum distributed load q_0: {q_0}")
print("\nThe force F is calculated such that the total deflection at the tip is zero.")
print("The deflection from the distributed load q(x) cancels the deflection from the point force F.")
print(f"Equation for F: F = (13 * q_0 * L) / 270")
print(f"Substituting values: F = (13 * {q_0} * {L}) / 270")
print(f"Final calculation: F = {13*q_0*L} / 270 = {F}")
print("\nThe required force F is:")
print(F)
print("<<<13.0>>>")
