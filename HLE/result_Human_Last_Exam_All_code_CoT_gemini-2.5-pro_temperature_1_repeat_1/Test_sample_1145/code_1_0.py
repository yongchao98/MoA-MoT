import numpy as np
from scipy.integrate import quad

# 1. Define the physical parameters for the problem.
#    - n: the quantum number
#    - a: the length of the box (we can set a=1.0 as the probability is dimensionless)
#    - x1, x2: the boundaries of the interval
n = 2
a = 1.0
x1 = 0.495 * a
x2 = 0.505 * a

# 2. Define the probability density function |ψ(x)|², which is (2/a) * sin²(n*π*x/a).
def probability_density(x, n_val, a_val):
    """
    Calculates the probability density |ψ(x)|² for a particle in a 1D box.
    """
    return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# 3. Use the `scipy.integrate.quad` function to calculate the definite integral.
#    This integral gives the probability. The function returns the result and an error estimate.
probability, error_estimate = quad(probability_density, x1, x2, args=(n, a))

# 4. Print the parameters, the symbolic representation of the integral, and the final result.
print("--- Calculation Details ---")
print(f"Quantum state n = {n}")
print(f"Box length a = {a}")
print(f"Interval = from {x1}a to {x2}a")
print("-" * 27)
print("The equation being solved is:")
print(f"P = Integral from {x1} to {x2} of (2/{a}) * sin^2({n}*pi*x/{a}) dx")
print("-" * 27)
print(f"The calculated probability is: {probability}")