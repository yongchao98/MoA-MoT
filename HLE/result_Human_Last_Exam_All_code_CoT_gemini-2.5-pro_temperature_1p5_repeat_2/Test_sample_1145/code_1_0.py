import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    within a specific interval for n=2.
    """
    # 1. Define the parameters for the problem.
    # The quantum state is n.
    n = 2
    # The interval is from 0.495a to 0.505a.
    # For calculation, we can set the box length 'a' to 1, as it cancels out.
    a = 1.0
    x1 = 0.495 * a
    x2 = 0.505 * a

    # 2. Define the probability density function, |ψ(x)|².
    def probability_density(x, n_val, a_val):
        """
        The probability density for a particle in a 1D box: (2/a) * sin²(nπx/a).
        """
        return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # 3. Integrate the probability density over the interval [x1, x2] to find the probability.
    # The quad function returns the integral result and an error estimate.
    probability, error = quad(probability_density, x1, x2, args=(n, a))

    # 4. Print the final equation and the result.
    # The prompt requires showing the numbers in the final equation.
    print("The probability (P) is calculated by integrating the probability density function |ψ(x)|².")
    print("The equation with the given values is:")
    print(f"P = Integral from {x1:.3f}a to {x2:.3f}a of (2/a) * sin^2({n}*pi*x/a) dx")
    print("\nAfter performing the integration, the result is:")
    print(f"Probability = {probability:.9f}")

solve_probability()