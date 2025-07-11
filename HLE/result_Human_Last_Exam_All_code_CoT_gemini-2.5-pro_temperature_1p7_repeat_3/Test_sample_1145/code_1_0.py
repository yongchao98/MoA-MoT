import numpy as np
from scipy.integrate import quad

# Define constants from the problem
n = 2
x1_factor = 0.495
x2_factor = 0.505
# We can set the box length 'a' to 1.0 for simplicity. 
# The probability is a dimensionless ratio, so the result is independent of 'a'.
a = 1.0

# Explain the problem and the physics
print("To find the probability of the particle in a 1D box between two positions,")
print("we integrate the probability density function |ψ(x)|² over that interval.")
print("The probability density is |ψ_n(x)|² = (2/a) * sin²(n*pi*x/a).\n")

# State the specific equation being solved, showing the numbers.
print("For this problem, the constants are:")
print(f"Quantum state, n = {n}")
print(f"Interval = [{x1_factor}a, {x2_factor}a]\n")

print("The equation for the probability (P) is the integral:")
print(f"P = Integral from x={x1_factor}a to x={x2_factor}a of (2/a) * sin^2({n}*pi*x/a) dx\n")

# Define the function to integrate
def probability_density(x, n_val, L):
    """
    Calculates the probability density |ψ(x)|² for a particle in a 1D box.
    """
    return (2.0 / L) * (np.sin(n_val * np.pi * x / L))**2

# Calculate the integration limits for a=1
x1 = x1_factor * a
x2 = x2_factor * a

# Perform the numerical integration
# The args=(n, a) passes the additional parameters (n_val, L) to the function.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# Print the final answer
print(f"The calculated probability is:")
print(probability)
print(f"(The estimated numerical error is: {error})")