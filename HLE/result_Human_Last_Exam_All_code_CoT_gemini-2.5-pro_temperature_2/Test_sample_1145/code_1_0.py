import numpy as np
from scipy import integrate

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box in a given interval.
    """
    # Define the parameters for the problem
    n = 2  # Quantum number
    # We can set the box length 'a' to 1.0, as it will cancel out mathematically.
    a = 1.0
    x1 = 0.495 * a
    x2 = 0.505 * a

    # Define the probability density function, |Ψ(x)|², for the particle in a 1D box.
    def probability_density(x, n_val, a_val):
        """
        The function to be integrated: (2/a) * sin^2(n*pi*x/a).
        """
        return (2 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # Perform the numerical integration using scipy.integrate.quad
    # quad returns the integral result and an estimated error.
    probability, _ = integrate.quad(
        probability_density,
        x1,
        x2,
        args=(n, a)
    )

    # Print the equation being solved, with all numbers included
    print("The probability (P) is found by solving the following integral:")
    print(f"P = integral from {x1:.3f} to {x2:.3f} of (2/{a:.1f}) * sin^2({n}*pi*x/{a:.1f}) dx")
    
    # Print the final result
    print(f"\nThe calculated probability is: {probability}")

if __name__ == "__main__":
    solve_probability()
