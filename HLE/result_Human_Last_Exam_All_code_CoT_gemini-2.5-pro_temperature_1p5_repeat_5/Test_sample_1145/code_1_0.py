import numpy as np
from scipy.integrate import quad

# This script calculates the probability of finding a particle in a 1D box
# within a given interval. If you don't have scipy installed, please run:
# pip install scipy

# 1. Define the parameters for the problem
n = 2           # The principal quantum number
a = 1.0         # We can set the box length 'a' to 1.0, as it cancels out
x1 = 0.495 * a  # The lower bound of the interval
x2 = 0.505 * a  # The upper bound of the interval

# 2. Define the probability density function, |ψ(x)|²
def probability_density(x, n, a):
  """
  Calculates the probability density |ψ(x)|² for a particle in a 1D box.
  """
  return (2/a) * (np.sin(n * np.pi * x / a))**2

# 3. Perform the numerical integration
#    quad returns a tuple: (result, estimated_error)
probability, error = quad(probability_density, x1, x2, args=(n, a))

# 4. Print the final result and the equation used
print("The probability P is found by integrating the probability density |ψ_n(x)|² over the interval [x1, x2]:")
print("P = Integral from x1 to x2 of (2/a) * sin²(n * π * x / a) dx")
print("\nFor this problem, the values are:")
print(f"n = {n}")
print(f"a = {a}")
print(f"x1 = {x1}")
print(f"x2 = {x2}")
print("\nThe final calculated probability is:")
print(f"P = {probability:.4e}")

# The problem is centered around a node (at x=a/2 for n=2), so the probability is very small.

print(f"<<<{probability}>>>")