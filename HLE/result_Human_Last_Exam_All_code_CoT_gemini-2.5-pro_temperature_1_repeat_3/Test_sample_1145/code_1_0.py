import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    in a given interval for a specific quantum state n.
    """
    # Define the parameters for the calculation
    n = 2
    x1_norm = 0.495  # Lower bound of the interval (x/a)
    x2_norm = 0.505  # Upper bound of the interval (x/a)

    # The equation to solve is the integral of the probability density function.
    # The probability density as a function of u = x/a is 2 * sin^2(n*pi*u).
    print(f"The probability P is the integral of 2 * sin^2({n}*pi*u) du from u={x1_norm} to u={x2_norm}.")

    # Define the integrand function
    def probability_density(u):
        return 2 * (np.sin(n * np.pi * u))**2

    # Perform the numerical integration using scipy.integrate.quad
    # quad returns a tuple: (result, estimated_error)
    probability, error = quad(probability_density, x1_norm, x2_norm)

    # Print the final calculated probability
    print(f"The calculated probability is: {probability}")

solve_probability()