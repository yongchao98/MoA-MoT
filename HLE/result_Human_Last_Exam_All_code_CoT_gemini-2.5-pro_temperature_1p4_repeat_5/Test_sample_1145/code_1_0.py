import numpy as np
from scipy.integrate import quad

# 1. Define the constants for the problem.
n = 2
x1_factor = 0.495
x2_factor = 0.505

# For the numerical calculation, we can set the length of the box 'a' to 1.
# The result is a dimensionless probability, as the 'a' dependence cancels out.
a = 1.0
x1 = x1_factor * a
x2 = x2_factor * a

# 2. Define the probability density function |ψ(x)|² as a Python function.
# The function is (2/a) * sin²(n*pi*x/a).
def probability_density(x, n_val, a_val):
  """
  Calculates the probability density for a particle in a 1D box.
  """
  return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# 3. Perform the numerical integration using scipy.integrate.quad.
# The quad function returns the integral result and an estimated error. We only need the result.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# 4. Print the final result, showing the equation being solved as requested.
print(f"The probability (P) is found by integrating the probability density function |ψ_n(x)|².")
print(f"For n={n}, the integral to solve is:")
print(f"P = Integral from {x1_factor}a to {x2_factor}a of (2/a) * sin²({n}πx/a) dx")
print("\nAfter performing the integration, the result is:")
print(f"Probability = {probability}")
