import numpy as np
from scipy.integrate import quad

# Define the parameters for the particle in a 1D box problem
n = 2  # The principal quantum number
a = 1.0  # We can set the length of the box 'a' to 1 for simplicity
x1 = 0.495 * a  # The start of the interval
x2 = 0.505 * a  # The end of the interval

def probability_density(x, n_val, a_val):
    """
    This function calculates the probability density |ψ(x)|² for a particle
    in a 1D box of length 'a' and quantum state 'n'.
    The formula is (2/a) * sin^2(n*pi*x/a).
    """
    return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# Perform the numerical integration of the probability density function
# over the interval [x1, x2].
# quad returns the result of the integral and an error estimate.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# The problem asks to output the numbers in the final equation.
# We will print a clear statement showing the parameters used for the calculation
# and the resulting probability.
print(f"The equation to solve is the integral of (2/a) * sin^2(n*pi*x/a) dx.")
print(f"For quantum number n = {n}, and box length a:")
print(f"The probability of finding the particle between x = {x1/a}a and x = {x2/a}a is:")
print(f"P = {probability}")