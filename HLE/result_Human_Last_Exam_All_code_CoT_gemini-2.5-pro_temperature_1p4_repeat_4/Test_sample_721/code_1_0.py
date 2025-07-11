import numpy as np
from scipy.integrate import quad
import math

# For function 2, we test the inequality sum(n*|a_n|^2) <= sum(|a_n|)
# by computing the LHS (which is Area/pi) and a lower bound for the RHS.

print("--- Analysis for Function 2 ---")
# The image of the unit disk under f(z) is a square. We first compute its side length K.
def integrand_K(x):
    """Integrand to compute the side length K of the square image."""
    # This integral is improper at x=0 and x=1 but converges.
    return 1.0 / np.sqrt(x * (1.0 - x**2))

# Use scipy's quad to perform the numerical integration.
K, K_err = quad(integrand_K, 0, 1)

# The area of the square image is K^2.
Area = K**2

# The LHS of the inequality is Area/pi.
LHS = Area / np.pi

# Now we compute a lower bound for the RHS. RHS = sum(|a_n|).
# RHS >= |a_0| + |a_1|.

# |a_0| is given by another integral.
def integrand_a0(t):
    """Integrand to compute |a_0| = |f(0)|."""
    return 1.0 / np.sqrt(t * (1.0 + t**2))

abs_a0, a0_err = quad(integrand_a0, 0, 1)

# |a_1| = |f'(0)| can be calculated analytically to be sqrt(2).
abs_a1 = np.sqrt(2)

# Lower bound for the RHS.
RHS_lower_bound = abs_a0 + abs_a1

# The final inequality check is LHS <= RHS. We check if LHS <= RHS_lower_bound.
print("The inequality to check is: sum(n*|a_n|^2)  <=  sum(|a_n|)")
print("This is equivalent to:       Area(f(D))/pi  <=  sum(|a_n|)\n")

print(f"The left-hand side (LHS) is Area/pi.")
print(f"Numerically, LHS = {LHS:.4f}")

print("\nThe right-hand side (RHS) is sum(|a_n|).")
print(f"A lower bound for the RHS is |a_0| + |a_1|.")
print(f"|a_0| = {abs_a0:.4f}")
print(f"|a_1| = {abs_a1:.4f}")
print(f"So, RHS >= {abs_a0:.4f} + {abs_a1:.4f} = {RHS_lower_bound:.4f}")

print("\n--- Final Check ---")
print(f"We check if the inequality {LHS:.4f} <= {RHS_lower_bound:.4f} holds.")
is_satisfied = LHS <= RHS_lower_bound
print(f"Result: {is_satisfied}")
print("\nSince the LHS is less than or equal to a lower bound of the RHS, the inequality holds for Function 2.")
