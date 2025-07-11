import numpy as np
from scipy.integrate import quad

# 1. Define the parameters for the problem based on the user's request.
n = 2  # The quantum number for the energy state.
a = 1.0  # We can set the box length 'a' to 1.0 for simplicity, as it cancels out.
x1_ratio = 0.495 # The starting point of the interval as a fraction of 'a'.
x2_ratio = 0.505 # The ending point of the interval as a fraction of 'a'.
x1 = x1_ratio * a # The lower bound of the interval.
x2 = x2_ratio * a # The upper bound of the interval.

# 2. Define the probability density function, P(x) = |ψ(x)|².
#    ψ_n(x) = sqrt(2/a) * sin(n*pi*x/a)
#    P(x) = (2/a) * sin^2(n*pi*x/a)
def probability_density(x, n_val, a_val):
    """
    This function represents the probability density for a particle in a 1D box.
    """
    return (2 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

# 3. Use numerical integration (scipy.integrate.quad) to find the probability.
#    The quad function returns the integral result and an estimated error.
#    We only need the first value, which is the result of the integration.
probability, error = quad(probability_density, x1, x2, args=(n, a))

# 4. Print the final result, showing the parameters used in the calculation as requested.
print(f"The equation to solve is the integral of the probability density function:")
print(f"P = Integral from {x1_ratio}a to {x2_ratio}a of (2/a) * sin^2({n}*pi*x/a) dx")
print("\n--- Calculation Result ---")
print(f"The probability of finding the particle in the state n = {n} between x = {x1_ratio}a and x = {x2_ratio}a is:")
print(f"{probability}")