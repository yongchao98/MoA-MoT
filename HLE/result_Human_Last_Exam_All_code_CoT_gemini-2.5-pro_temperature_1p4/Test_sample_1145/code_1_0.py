import numpy as np
from scipy.integrate import quad

# Define the constants for the problem.
# n is the principal quantum number.
n = 2
# We can set the length of the box 'a' to 1.0, as the 'a' units cancel out
# in the final probability calculation, making the result dimensionless.
a = 1.0
# Define the interval [x1, x2] as a fraction of the box length 'a'.
x1 = 0.495 * a
x2 = 0.505 * a

# The probability density function is the square of the normalized wave function.
# For a particle in a 1D box: |Ψ(x)|² = (2/a) * sin²(n*π*x/a).
def probability_density(x, n_val, a_val):
    """Calculates the probability density for a particle in a 1D box."""
    return (2 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# To find the probability of the particle being in the interval [x1, x2],
# we integrate the probability density function over this interval.
# The quad function returns the integral result and an estimate of the error.
probability, _ = quad(probability_density, x1, x2, args=(n, a))

# --- Output the final equation and the result ---

# Print the parameters of the final equation
print(f"The equation for the probability (P) is the integral of the probability density function |Ψ(x)|².")
print(f"For n = {n}, the function is (2/a) * sin²({n}*π*x/a).")
print(f"The integral is evaluated from x₁ = {x1:.3f}a to x₂ = {x2:.3f}a.")
print("\nFinal Equation:")
print(f"P = ∫ from {x1:.3f} to {x2:.3f} of (2/{a:.1f}) * sin²({n}*π*x/{a:.1f}) dx")

# Print the final calculated probability
print(f"\nThe calculated probability is: {probability}")