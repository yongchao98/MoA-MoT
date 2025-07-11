import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box in a given interval.
    """
    # Define the parameters of the problem
    n = 2
    a = 1.0  # The length of the box 'a' is set to 1 for calculation.
    x1 = 0.495 * a
    x2 = 0.505 * a

    # The probability density function for a particle in a 1D box
    # |psi_n(x)|^2 = (2/a) * sin^2(n*pi*x/a)
    def probability_density(x, n_val, L):
        return (2.0 / L) * (np.sin(n_val * np.pi * x / L))**2

    # Calculate the probability by integrating the probability density function
    # over the specified interval [x1, x2].
    probability, _ = quad(probability_density, x1, x2, args=(n, a))

    # --- Output Section ---
    print("The probability P is found by integrating the probability density function |psi_n(x)|^2 from x1 to x2.")
    
    # Print the equation with all the specific numbers, as requested.
    # We use a placeholder 'a' for clarity, then show the version with a=1 used for the calculation.
    print("\nThe equation to solve for this problem is:")
    print(f"P = integral from x1={0.495}a to x2={0.505}a of (2/a) * sin^2({n}*pi*x/a) dx")

    print("\nFor the numerical calculation, we can set a = 1.0:")
    print(f"P = integral from {x1} to {x2} of ({2.0}/{a}) * sin^2({n}*pi*x/{a}) dx")
    
    print(f"\nThe calculated probability is: {probability:.10f}")

solve_probability()