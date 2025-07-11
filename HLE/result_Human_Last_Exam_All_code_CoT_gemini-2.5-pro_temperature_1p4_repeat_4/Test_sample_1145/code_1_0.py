import numpy as np
from scipy.integrate import quad

# --- Problem Parameters ---
# Quantum number
n = 2
# Length of the box (we can set a=1 as it cancels out)
a = 1.0
# Interval limits
x1 = 0.495 * a
x2 = 0.505 * a

# The probability density function is |Ψ(x)|², where Ψ(x) is the wavefunction.
# For a particle in a 1D box: |Ψ(x)|² = (2/a) * sin²(n*π*x/a)
def probability_density(x, n_val, a_val):
    """Calculates the probability density for a particle in a 1D box."""
    return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# Calculate the probability by integrating the density function over the interval [x1, x2].
# The quad function from the scipy library performs numerical integration.
# It returns the result of the integral and an estimate of the error.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# --- Output the results ---
# As requested, we will print the final equation with all numbers included.
print("To find the probability, we integrate the probability density function P(x) = |Ψ(x)|².")
print("The equation for the probability is:")
# The formatting below shows the integral with the specific values for this problem.
print(f"P = ∫[from {x1:.3f}a to {x2:.3f}a] (2/a) * sin²({n}πx/a) dx")
print("\nCalculating this integral gives the probability.")
print(f"The probability is: {probability}")