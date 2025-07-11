import numpy as np
from scipy.integrate import quad

# --- Step 1: Define given constants and calculate geometric properties ---
# The parameter 'a' is given
a = 12**(1/4)
a_4 = a**4 # This will be 12.0

# Calculate the moments of inertia I_ss and I_zz for the given cross-section.
# The cross-section is a large square (3a x 3a) with two smaller squares (a x a) removed.
# Using the parallel axis theorem, the moments of inertia for the final shape are:
# I_zz = I_zz_large - 2 * I_zz_cutout
# I_ss = I_ss_large - 2 * I_ss_cutout
# After calculation, the simplified formulas are:
I_zz = (73 * a_4) / 12
I_ss = (55 * a_4) / 12

# --- Step 2: Calculate the problem parameters L, q0, and EI ---
# These are given in terms of I_zz, I_ss, and an integral.
L = (30 * I_zz) / 73
q0 = (9 * I_ss) / 55

# The bending stiffness EI is given by the Fresnel sine integral.
def integrand(x):
    return np.sin(x**2)
EI, _ = quad(integrand, 0, np.pi)

# --- Step 3: Formulate and solve for the force F ---
# From the principle of superposition, the total deflection at the tip is zero.
# y_total = y_F + y_q = 0
# (9 * F * L^3 / (8 * EI)) - (13 * q0 * L^4 / (240 * EI)) = 0
# Solving for F gives the formula:
# F = (13 * q0 * L) / 270

# --- Step 4: Calculate the final numerical value of F and print the results ---
F_final = (13 * q0 * L) / 270

print("Step 1: Calculate geometric and material properties")
print(f"a = {a:.4f}")
print(f"I_zz = {I_zz:.4f}")
print(f"I_ss = {I_ss:.4f}")
print(f"L = {L:.4f}")
print(f"q_0 = {q0:.4f}")
print(f"EI = {EI:.4f}")
print("-" * 30)
print("Step 2: Solve for the force F")
print("The equation for the force F is derived from setting the total deflection at the tip to zero:")
print("F = (13 * q_0 * L) / 270")
print("\nSubstituting the calculated values:")
# The final equation with each number
print(f"F = (13 * {q0:.2f} * {L:.2f}) / 270")
print(f"F = {F_final:.4f}")

print("\n<<<13.0>>>")