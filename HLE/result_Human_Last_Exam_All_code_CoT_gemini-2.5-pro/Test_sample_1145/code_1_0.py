import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    in a given interval for a specific quantum number n.
    """
    # Define the parameters for the calculation
    # n is the principal quantum number.
    n = 2
    
    # The integration limits are given as 0.495a and 0.505a.
    # By using a change of variables u = x/a, the box length 'a' cancels out,
    # so we can integrate from u=0.495 to u=0.505.
    x1 = 0.495
    x2 = 0.505

    # The probability density function is |ψ(u)|² = 2 * sin²(n*π*u), where u=x/a.
    def probability_density(u, n_level):
        return 2 * (np.sin(n_level * np.pi * u))**2

    # Calculate the definite integral of the probability density function
    # over the interval [x1, x2].
    # The quad function returns the integration result and an error estimate.
    probability, error = quad(probability_density, x1, x2, args=(n,))

    # Print the equation and the values used for the calculation,
    # followed by the final result.
    print("The probability P is found by integrating the probability density |ψ_n(x)|² from x1 to x2.")
    print("The final equation is: P = integral from x1 to x2 of (2/a) * sin^2(n*π*x/a) dx")
    print("\nFor the given values:")
    print(f"Quantum number n = {n}")
    print(f"Lower bound x1 = {x1}a")
    print(f"Upper bound x2 = {x2}a")
    print("\nThe calculated probability is:")
    print(f"P = {probability}")

solve_probability()
<<<6.579708171452109e-06>>>