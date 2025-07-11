import numpy as np
from scipy.integrate import quad

# The final equation to integrate is of the form:
# Integral from 0 to infinity of a / (b + c * exp(d * tau^2)) d_tau
# We determined the coefficients to be:
a = 4.0
b = 1.0
c = 3.0
d = 2.0

# Print the numbers in the final equation
print(f"The integral to be computed is the integral from 0 to infinity of:")
print(f"{a} / ({b} + {c} * exp({d} * tau^2)) d_tau")
print("-" * 20)
print(f"Coefficients:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"d = {d}")
print("-" * 20)


def integrand(tau):
  """The function to be integrated."""
  return a / (b + c * np.exp(d * tau**2))

# Perform the numerical integration from 0 to infinity.
result, error = quad(integrand, 0, np.inf)

# Print the final result
print(f"The result of the integration is: {result}")
print(f"The estimated error is: {error}")

# The final answer in the requested format
print(f"\n<<<{result}>>>")
