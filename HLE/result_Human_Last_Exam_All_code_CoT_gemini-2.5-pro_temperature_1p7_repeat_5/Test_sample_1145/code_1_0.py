import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    within a given interval for a specific quantum state.
    """
    # Define the constants for the problem based on the user's request.
    # Quantum state level
    n = 2
    # The interval for x is from 0.495a to 0.505a.
    # We use a change of variables u = x/a, so the integration limits for u are:
    u1 = 0.495
    u2 = 0.505

    # The function to integrate is the probability density after the change of variables.
    # Integrand = 2 * sin^2(n * pi * u)
    def integrand(u, n_val):
        return 2 * (np.sin(n_val * np.pi * u))**2

    # Perform the numerical integration using scipy.integrate.quad.
    # The quad function returns the integral result and an estimate of the error.
    probability, error = quad(integrand, u1, u2, args=(n,))

    # Output the parameters of the problem and the final result.
    # The user requested to output each number in the final equation.
    print(f"Problem: Find the probability of a particle in a 1D box.")
    print(f"The equation for probability is the integral of 2 * sin^2(n * pi * x / a) dx / a")
    print(f"For the quantum state n = {n},")
    print(f"in the interval from x = {u1}*a to {u2}*a.")
    
    # Print the final calculated probability.
    print("\nResult:")
    print(f"The probability is: {probability}")

if __name__ == "__main__":
    solve_probability()
