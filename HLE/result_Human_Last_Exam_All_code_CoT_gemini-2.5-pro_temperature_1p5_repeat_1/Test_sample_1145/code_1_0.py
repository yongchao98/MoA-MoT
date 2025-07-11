import numpy as np
from scipy.integrate import quad

# Define the constants for the problem
n = 2
lower_bound_ratio = 0.495
upper_bound_ratio = 0.505

# The function to integrate is the probability density after a change of variables u = x/a.
# P(u)du = 2 * sin²(n*pi*u) du
def probability_density(u, n_val):
    """
    Calculates the probability density as a function of the normalized position u=x/a.
    """
    return 2 * (np.sin(n_val * np.pi * u))**2

# Perform the numerical integration using scipy.integrate.quad.
# The quad function returns the integral result and an estimated error.
probability, error = quad(
    probability_density,
    lower_bound_ratio,
    upper_bound_ratio,
    args=(n,) # Pass the value of n to the function
)

# Output the explanation, the final equation with all numbers, and the result.
print(f"To find the probability of the particle's location, we integrate the probability density function.")
print(f"For n = {n}, from x = {lower_bound_ratio}a to {upper_bound_ratio}a, the equation to solve is:")
print(f"Probability = ∫ from {lower_bound_ratio}a to {upper_bound_ratio}a of (2/a) * sin²({n}πx/a) dx")
print(f"\nThis integral is calculated numerically. The resulting probability is:")
print(f"{probability}")