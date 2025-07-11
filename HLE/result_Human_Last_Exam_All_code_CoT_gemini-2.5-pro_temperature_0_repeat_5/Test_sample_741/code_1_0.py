import scipy.special as sc
from scipy.optimize import fsolve

# The summation is equal to the modified Bessel function I_{x-1}(2).
# We need to find the largest x for which I_{x-1}(2) = 0.
# Let nu = x - 1. We need to find the largest root of I_nu(2) = 0.

# Define the function for the root-finding algorithm.
def bessel_function(nu):
  """
  Computes the modified Bessel function of the first kind, I_nu(2).
  """
  return sc.iv(nu, 2)

# From mathematical properties and plots, the largest root is negative,
# located near -2.3. We use this as an initial guess for fsolve.
initial_guess_nu = -2.3
largest_root_nu = fsolve(bessel_function, initial_guess_nu)[0]

# Calculate x from the relationship nu = x - 1
x = largest_root_nu + 1

# Output the steps of the calculation
print(f"The equation to solve is I_(x-1)(2) = 0.")
print(f"Let nu = x - 1. We solve for the largest root of I_nu(2) = 0.")
print(f"The largest root found is nu = {largest_root_nu}")
print(f"The value of x is calculated as x = nu + 1.")
print(f"x = {largest_root_nu} + 1 = {x}")
print(f"The final answer in the format {{-a.bbb}} is:")
print(f"{{{x:.3f}}}")
