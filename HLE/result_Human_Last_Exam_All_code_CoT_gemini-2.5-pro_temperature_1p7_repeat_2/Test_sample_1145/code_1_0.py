import numpy as np
from scipy.integrate import quad

# 1. Define the constants for the problem.
# We set the box length 'a' to 1 because the final probability is independent of it.
n = 2
a = 1.0
x1 = 0.495 * a
x2 = 0.505 * a

# 2. Define the probability density function to be integrated.
# This is the square of the magnitude of the wavefunction: |ψ_n(x)|²
def probability_density(x, n_val, a_val):
    """
    Probability density function for a particle in a 1D box.
    P(x) = (2/a) * sin^2(n*pi*x/a)
    """
    return (2 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# 3. Calculate the probability by integrating the density function over the interval [x1, x2].
# The quad function returns the integral result and an estimated error, we only need the result.
probability, _ = quad(probability_density, x1, x2, args=(n, a))

# 4. Print the final equation with all the numbers and the result.
# The equation is shown for a=1, as used in the calculation.
print("The probability is calculated by integrating the probability density function P(x) over the interval.")
print("The final equation with the specified values (using a=1 for calculation) is:")
print(f"Probability = Integral from {x1} to {x2} of (2/{a}) * sin^2({n} * pi * x / {a}) dx")
print("\nResult:")
print(f"The probability that the particle is located in the interval from {x1}a to {x2}a for n={n} is: {probability}")
print(f"In scientific notation, the probability is approximately: {probability:.3e}")
