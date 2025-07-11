import scipy.special as sp
import scipy.optimize as opt
import numpy as np

# The problem is to find the largest x for which the summation equals 0.
# The summation is equivalent to the modified Bessel function I_{x-1}(2).
# We need to solve the equation I_{x-1}(2) = 0.
# Let alpha = x - 1. We find the largest root 'alpha' of I_alpha(2) = 0.
# It is known that roots exist only for alpha < -1.
# By checking values, we find the largest root lies between -2.6 and -2.5.

def bessel_function_at_2(alpha):
  """Defines the function I_alpha(2) whose root we need to find."""
  return sp.iv(alpha, 2.0)

# Use a numerical root-finding algorithm to find the largest root for alpha.
try:
    solution = opt.root_scalar(bessel_function_at_2, bracket=[-2.6, -2.5], method='brentq')
    alpha_root = solution.root

    # Calculate the corresponding x value using the relation x = alpha + 1.
    x_value = alpha_root + 1

    # Print the final equation with the computed numbers, formatted to three decimal places.
    print("The largest root 'alpha' of the equation I_alpha(2) = 0 is found to be:")
    print(f"alpha = {alpha_root:.3f}")
    print("\nThe value of x is calculated from the relation x = alpha + 1:")
    print(f"x = {alpha_root:.3f} + 1 = {x_value:.3f}")
    print("\nThe largest x value, in the required format, is:")
    print(f"{{{x_value:.3f}}}")

except (ImportError, ModuleNotFoundError):
    print("Please install scipy and numpy libraries to run this code.")
    print("You can install them using: pip install scipy numpy")
