import numpy as np
from scipy.integrate import quad

# --- Problem Parameters ---
# Quantum number
n = 2
# Length of the box 'a'. It cancels out, so we set it to 1.0 for calculation.
a = 1.0
# Interval boundaries
x1 = 0.495 * a
x2 = 0.505 * a

# --- Probability Density Function ---
# The function to be integrated: |psi_n(x)|^2 = (2/a) * sin^2(n*pi*x/a)
def probability_density(x, n_val, a_val):
    """Calculates the probability density for a particle in a 1D box."""
    return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# --- Calculation ---
# Use scipy.integrate.quad to perform the numerical integration.
# It returns the integral result and an estimate of the absolute error.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# --- Output ---
# Print the final equation showing all the numbers and the calculated result.
# The result is formatted in scientific notation for clarity.
print("The probability is calculated by the integral:")
print(f"P = Integral from {x1} to {x2} of (2/{a})*sin^2({n}*pi*x/{a}) dx = {probability:.6e}")