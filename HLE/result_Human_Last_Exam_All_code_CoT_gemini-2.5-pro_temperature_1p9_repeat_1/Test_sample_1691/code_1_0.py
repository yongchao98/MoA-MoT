import numpy as np

# This script calculates the parameters for the analytical approximation
# of the integral I(epsilon) for small epsilon.
# The approximation has the form: I(epsilon) ≈ C * epsilon^p.

# From the problem, the integrand is 1/(epsilon + 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0).
# For small x, the dominant term in the denominator is a*x^n.
a = 9.0
n = 5.0

# ===== Calculation =====
# The exponent 'p' in the approximation formula is p = (1/n) - 1.
p_val = (1.0 / n) - 1.0

# The coefficient 'C' is given by C = (a**(-1/n)) * (pi / (n * sin(pi/n))).
C_val = (a**(-1.0 / n)) * (np.pi / (n * np.sin(np.pi / n)))


# ===== Output =====
# The prompt requests that we output each number in the final equation.
# The final equation is I(epsilon) ≈ C * epsilon^p

# We will print the equation with its numerical coefficients.
print("The analytical formula that approximates the integral is:")
print(f"I(epsilon) ≈ {C_val} * epsilon^({p_val})")
print("\nWhere each number in the final equation is:")
print(f"The coefficient C = {C_val}")
print(f"The exponent p = {p_val}")