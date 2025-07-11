import numpy as np
from scipy.integrate import quad

# Define the constants for the problem
n = 2
# We can set the length of the box 'a' to 1 for simplicity, as it cancels out.
a = 1.0
x1 = 0.495 * a
x2 = 0.505 * a

# Define the probability density function P(x) = |ψ(x)|^2
# P(x) = (2/a) * sin^2(n*π*x/a)
def probability_density(x, n_val, a_val):
    """
    Calculates the probability density for a particle in a 1D box.
    """
    return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# --- Calculation ---
# The probability is the integral of the probability density function
# over the specified interval.

# Print the equation being solved
print("We are calculating the probability by solving the definite integral:")
print(f"Probability = Integral from {x1:.3f}a to {x2:.3f}a of (2/a) * sin^2({n}*pi*x/a) dx")
print("-" * 20)

# Perform the numerical integration using scipy.integrate.quad
# The quad function returns the result of the integral and an estimate of the error.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# Print the final result
print(f"The calculated probability is: {probability}")
<<<6.539811528180436e-06>>>