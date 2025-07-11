import numpy as np
from scipy.integrate import quad

# Step 1: Define the constants for the problem.
# The quantum number n.
n = 2
# The length of the box 'a'. We can set it to 1.0 since the 'a'
# dependency cancels out and the probability is dimensionless.
a = 1.0
# The lower bound of the interval.
x1 = 0.495 * a
# The upper bound of the interval.
x2 = 0.505 * a

# Step 2: Define the probability density function, |psi(x)|^2.
def probability_density(x, n_val, a_val):
  """
  Calculates the probability density for a particle in a 1D box.
  |psi(x)|^2 = (2/a) * sin^2(n*pi*x/a)
  """
  return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# Step 3: Calculate the probability by numerically integrating the
# probability density function over the interval [x1, x2].
# The quad function returns the result of the integral and an estimate of the error.
probability, error_estimate = quad(probability_density, x1, x2, args=(n, a))

# Step 4: Print the inputs and the final calculated probability.
print("Calculating the probability that the particle in a 1D box is located in an interval from x1 to x2.")
print(f"Quantum state (n): {n}")
print(f"Interval start (x1): {x1}a")
print(f"Interval end (x2): {x2}a")
print("\nThe probability is the integral of (2/a) * sin^2(n*pi*x/a) from x1 to x2.")
print(f"\nCalculated Probability = {probability}")
