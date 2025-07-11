import numpy as np
from scipy.integrate import quad

# 1. Define the parameters from the problem statement.
n = 2
# The length of the box 'a' cancels out in the final probability calculation,
# so we can set a=1.0 for simplicity.
a = 1.0
x1 = 0.495 * a
x2 = 0.505 * a

# 2. Define the probability density function, |ψ(x)|².
def probability_density(x, n_val, a_val):
    """
    Probability density function |ψ_n(x)|² for a particle in a 1D box.
    """
    return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# 3. Calculate the probability by integrating the probability density function
# from x1 to x2 using scipy.integrate.quad.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# 4. Print the final equation with all the numbers and the result.
print("The probability P is found by integrating the probability density function |ψ_n(x)|² over the interval.")
print("The final equation to solve is:")
print(f"P = ∫ from x1={x1:.3f}a to x2={x2:.3f}a of (2/a) * sin²({n}*π*x/a) dx")
print("\nThe result of this integral is:")
print(f"{probability}")