import numpy as np
from scipy.integrate import quad

# For function 2, the inequality boils down to K <= pi / sqrt(2),
# where K is the value of a specific definite integral.

# Define the integrand for K
def integrand(t):
    return np.sin(t)**(-0.5)

# Calculate K by numerical integration from 0 to pi/2
# The integral has a singularity at 0, but it is integrable.
# We can compute it from a small epsilon to pi/2, or let quad handle it.
K, K_error = quad(integrand, 0, np.pi/2)

# Calculate pi / sqrt(2)
pi_over_sqrt2 = np.pi / np.sqrt(2)

# Calculate LHS and RHS of the original inequality
lhs = K**2 / np.pi
rhs = K / np.sqrt(2)

# The final inequality is K <= pi / sqrt(2)
# The numbers in the final equation are K, pi, and sqrt(2)
# Here we print the comparison between K and pi/sqrt(2)
print("Analysis for Function 2:")
print(f"The constant K is the integral of (sin(t))^(-1/2) from 0 to pi/2.")
print(f"Numerically calculated K = {K:.6f}")
print(f"The value to compare against is pi / sqrt(2) = {pi_over_sqrt2:.6f}")

print("\nThe inequality to check is K <= pi / sqrt(2).")
print(f"Check: {K:.6f} <= {pi_over_sqrt2:.6f}")
is_satisfied_2 = K <= pi_over_sqrt2
print(f"Is the inequality satisfied for Function 2? {is_satisfied_2}")

print("\nLet's check the original inequality form: sum(n*|a_n|^2) <= sum(|a_n|)")
print(f"LHS = K^2 / pi = {lhs:.6f}")
print(f"RHS = K / sqrt(2) = {rhs:.6f}")
print(f"Check: {lhs:.6f} <= {rhs:.6f}")
print(f"Is the inequality satisfied for Function 2? {lhs <= rhs}")

print("\nSummary of conclusions:")
print("Function 1: Does not satisfy the inequality (LHS diverges, RHS is finite).")
print("Function 2: Does not satisfy the inequality (as shown by numerical calculation).")
print("Function 3: Satisfies the inequality (LHS is finite, RHS diverges to infinity).")
print("\nTherefore, only function 3 satisfies the inequality.")