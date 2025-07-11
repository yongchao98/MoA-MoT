import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    in a given interval for a specific energy state.
    """
    # --- Parameters ---
    # The quantum number for the energy state
    n = 2
    
    # For calculation purposes, we can set the box length 'a' to 1.0.
    # The probability is a dimensionless quantity, and the 'a' dependency
    # cancels out when the interval is defined as a fraction of 'a'.
    a = 1.0
    
    # The interval factors are given in the problem
    x1_factor = 0.495
    x2_factor = 0.505
    
    # Calculate the absolute bounds of the integration
    x1 = x1_factor * a
    x2 = x2_factor * a

    # --- Probability Density Function ---
    # The probability density P(x) for a particle in a 1D box is |ψ(x)|²,
    # where ψ_n(x) = sqrt(2/a) * sin(n*π*x/a).
    # So, P(x) = (2/a) * sin²(n*π*x/a).
    def probability_density(x, n_val, a_val):
        return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # --- Calculation ---
    # To find the probability, we integrate the probability density function
    # over the interval [x1, x2] using scipy's quad function.
    probability, _ = quad(probability_density, x1, x2, args=(n, a))

    # --- Output ---
    # The prompt requires printing the numbers used in the final equation.
    print("The probability (P) of finding the particle in an interval is given by the integral of the probability density function:")
    print(f"P = Integral from x1 to x2 of (2/a) * sin²(nπx/a) dx")
    print("\nFor the specific problem, the values in the equation are:")
    print(f"Quantum number n = {n}")
    print(f"Lower bound x1 = {x1_factor}a")
    print(f"Upper bound x2 = {x2_factor}a")
    
    print(f"\nThe calculated probability is: {probability}")

solve_probability()