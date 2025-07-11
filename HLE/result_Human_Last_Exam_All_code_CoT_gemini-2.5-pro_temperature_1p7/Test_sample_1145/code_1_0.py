import numpy as np
from scipy.integrate import quad

# --- Parameters for the problem ---
# Quantum number
n = 2
# Length of the box (can be set to 1.0 as it cancels out)
a = 1.0
# Start of the interval
x1 = 0.495 * a
# End of the interval
x2 = 0.505 * a

# --- Probability Density Function ---
# P(x) = (2/a) * sin^2(n*pi*x / a)
# This is the function we will integrate.
def probability_density(x, n_val, a_val):
  """Calculates the probability density for a particle in a 1D box."""
  return (2 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# --- Numerical Integration ---
# Use scipy.integrate.quad to calculate the definite integral.
# quad returns two values: the result of the integration and an estimate of the error.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# --- Output the Final Result ---
# As requested, we print the final equation with the numbers and the result.
# The 'a' is left symbolically in the printout for clarity.
print("The probability P is found by solving the definite integral of the probability density function:")
print(f"P = Integral from {x1:.3f}a to {x2:.3f}a of (2/a) * sin^2({n}*pi*x/a) dx")
print("\nCalculating this integral gives the probability:")
print(f"P = {probability}")