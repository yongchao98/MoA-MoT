import numpy as np
from scipy.integrate import quad

def solve_particle_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    in a specific interval for a given energy state n.
    """
    # Define the constants for the problem
    n = 2
    x1_norm = 0.495
    x2_norm = 0.505
    
    # The length of the box 'a' can be set to 1.0 for simplicity,
    # as the final probability is independent of 'a'.
    a = 1.0
    
    # Define the integration limits
    x1 = x1_norm * a
    x2 = x2_norm * a

    # Define the probability density function, |ψ(x)|², to be integrated
    # |ψ_n(x)|² = (2/a) * sin²(n*π*x/a)
    def probability_density(x, n_val, a_val):
        return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # Calculate the definite integral using scipy.integrate.quad
    # The result is a tuple (value, estimated_error)
    probability, _ = quad(probability_density, x1, x2, args=(n, a))

    # Print the equation and the numbers involved
    print("The probability P is found by integrating the probability density |ψ_n(x)|² over the interval [x1, x2].")
    print(f"The equation for n = {n} from x1 = {x1_norm}a to x2 = {x2_norm}a is:")
    print(f"P = integral from {x1_norm}a to {x2_norm}a of (2/a) * sin^2({n}*pi*x/a) dx")
    print("\nCalculating this integral gives the probability:")
    print(f"P = {probability}")

# Execute the function
solve_particle_probability()
<<<6.579709664358983e-06>>>