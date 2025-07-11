import numpy as np
from scipy.integrate import quad

# The problem is to calculate the time-averaged integral:
# integral from 0 to infinity of 1 / (x(0; tau) + y(0; tau) + z(0; tau)) d(tau)
#
# Through analysis of the system's dynamics, we derived the expression for the sum
# of initial positions as a function of tau:
# x(0; tau) + y(0; tau) + z(0; tau) = (1/4) * (1 + 3 * exp(2 * tau^2))
#
# The integrand is the reciprocal of this expression:
# f(tau) = 4 / (1 + 3 * exp(2 * tau^2))

# Let's define the constants of the final equation for the integrand.
A = 4.0
B = 1.0
C = 3.0
D = 2.0

print("The integrand is of the form: A / (B + C * exp(D * tau**2))")
print(f"The constant values are:")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print(f"D = {D}")
print("-" * 20)

def integrand(tau):
  """
  The derived function to be integrated.
  """
  return A / (B + C * np.exp(D * tau**2))

# We use scipy.integrate.quad for numerical integration.
# The limits are from 0 to infinity (np.inf).
integral_value, error_estimate = quad(integrand, 0, np.inf)

print(f"The computed value of the integral is: {integral_value}")
print(f"The estimated absolute error is: {error_estimate}")

# The exact analytical result is (pi * sqrt(6)) / 9
# We can print this for comparison.
# analytical_result = (np.pi * np.sqrt(6)) / 9
# print(f"Analytical result for comparison: {analytical_result}")