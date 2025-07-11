import numpy as np
from scipy.integrate import quad

def probability_density(x, n, a):
    """
    Calculates the probability density |ψ(x)|² for a particle in a 1D box.
    
    Args:
        x (float): Position inside the box.
        n (int): The principal quantum number.
        a (float): The length of the box.
        
    Returns:
        float: The probability density at position x.
    """
    return (2 / a) * (np.sin(n * np.pi * x / a))**2

# --- Problem Parameters ---
# The quantum state
n = 2
# The interval of interest
x_start_factor = 0.495
x_end_factor = 0.505

# For calculation purposes, we can set the length of the box 'a' to 1.0.
# The probability is a dimensionless ratio, so 'a' will cancel out.
a = 1.0
x1 = x_start_factor * a
x2 = x_end_factor * a

# --- Calculation ---
# Use scipy.integrate.quad to perform the numerical integration.
# The 'args' tuple passes the parameters n and a to the probability_density function.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# --- Output ---
# Print the final equation with the specific numbers for this problem, as requested.
print("The probability P is found by solving the integral of the probability density |ψₙ(x)|² over the interval.")
print("The specific equation for this problem is:")
print(f"P = \u222B [from {x_start_factor}a to {x_end_factor}a] (2/a) * sin\u00b2({n}\u03c0x/a) dx\n")

print(f"The calculated probability for n={n} in the interval [{x1}, {x2}] is:")
print(probability)
