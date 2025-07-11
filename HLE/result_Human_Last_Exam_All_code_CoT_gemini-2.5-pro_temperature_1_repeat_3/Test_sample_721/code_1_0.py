import numpy as np
from scipy.integrate import quad
import math

# Plan:
# 1. For function 2, the map from the unit disk to a square, we analyze the inequality
#    Sum(n*|a_n|^2) <= Sum(|a_n|).
# 2. The LHS is Area/pi. The area is K1^2, where K1 is the side length of the square.
#    K1 = integral from 0 to 1 of 1/sqrt(x(1-x^2)) dx. With substitution u^2=x, this is
#    integral from 0 to 1 of 2/sqrt(1-u^4) du.
# 3. The RHS is Sum(|a_n|). A lower bound is |a0| + |a1|.
#    |a1| is sqrt(2).
#    |a0| = |integral from 0 to i of 1/sqrt(x(1-x^2)) dx|. With substitution x=it, this is
#    integral from 0 to 1 of 1/sqrt(t(1+t^2)) dt. With substitution u^2=t, this is
#    integral from 0 to 1 of 2/sqrt(1+u^4) du.
# 4. We will compute these values and compare LHS with the lower bound of RHS.

print("Verifying the inequality for Function 2: the conformal map to a square.")

def integrand_k1(u):
    """Integrand for the side length K1 of the square."""
    return 2 / np.sqrt(1 - u**4)

def integrand_abs_a0(u):
    """Integrand for the coefficient |a0|."""
    return 2 / np.sqrt(1 + u**4)

# Compute K1 and |a0|
K1, _ = quad(integrand_k1, 0, 1)
abs_a0, _ = quad(integrand_abs_a0, 0, 1)

# |a1| is known to be sqrt(2)
abs_a1 = np.sqrt(2)

# Calculate LHS and a lower bound for RHS
LHS = K1**2 / np.pi
RHS_lower_bound = abs_a0 + abs_a1

print("\n--- Calculation Results ---")
print(f"Side length of the square, K1 = {K1:.5f}")
print(f"Coefficient |a0| = {abs_a0:.5f}")
print(f"Coefficient |a1| = {abs_a1:.5f}")

print("\n--- Inequality Check ---")
print("The inequality is: Sum(n*|a_n|^2)  <=  Sum(|a_n|)")
print(f"LHS = Area/pi = K1^2/pi = ({K1:.5f})^2 / {np.pi:.5f} = {LHS:.5f}")
print(f"RHS >= |a0| + |a1| = {abs_a0:.5f} + {abs_a1:.5f} = {RHS_lower_bound:.5f}")

print("\n--- Comparison ---")
print(f"Is {LHS:.5f} <= {RHS_lower_bound:.5f}?  {LHS <= RHS_lower_bound}")
print("\nThe comparison shows that the LHS is smaller than a lower bound for the RHS.")
print("Therefore, the inequality holds for function 2.")