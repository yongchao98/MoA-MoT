import numpy as np
from scipy.integrate import quad

# Define the parameters for the calculation
n = 2
x1_factor = 0.495
x2_factor = 0.505

# The probability density function for a particle in a 1D box is P(x) = (2/a)*sin^2(n*pi*x/a).
# To calculate the probability of finding the particle in the interval [x1, x2],
# we integrate P(x) from x1 to x2.
# By substituting y = x/a, the integral becomes independent of 'a':
# Probability = Integral from y1 to y2 of 2*sin^2(n*pi*y) dy
# where y1 = x1/a and y2 = x2/a.

def integrand(y, n_val):
  """
  The integrand for the probability calculation, 2 * sin^2(n*pi*y).
  y is the normalized position x/a.
  n_val is the quantum state n.
  """
  return 2 * (np.sin(n_val * np.pi * y))**2

# Perform the numerical integration
# quad returns a tuple (integral_value, error_estimate)
probability, error = quad(integrand, x1_factor, x2_factor, args=(n,))

# Output the results, showing the equation being solved
print("For a particle in a 1D box of length 'a':")
print(f"The probability for the state n={n} in the interval [{x1_factor}a, {x2_factor}a] is calculated by the definite integral:")
print(f"P = Integral from {x1_factor}a to {x2_factor}a of (2/a) * sin^2({n}*pi*x/a) dx")
print("\nThis evaluates to:")
print(f"{probability}")