import numpy as np
from scipy.integrate import quad

def main():
    """
    Calculates the probability of finding a particle in a 1D box
    in a specific interval for a given quantum state n.
    """
    # Define the constants for the problem
    n = 2
    x1 = 0.495  # Lower bound of the interval (relative to box length a)
    x2 = 0.505  # Upper bound of the interval (relative to box length a)

    # The probability density function, |ψ_n(x)|², after setting a=1.
    # The 'a' term cancels out during integration.
    # integrand = 2 * sin²(nπx)
    def probability_density(x, n_val):
        return 2 * (np.sin(n_val * np.pi * x))**2

    # Perform the numerical integration
    # quad returns the result and an estimate of the absolute error
    probability, error = quad(probability_density, x1, x2, args=(n,))

    # The problem asks to output the numbers in the final equation.
    # The final equation we are solving is the simplified integral.
    print(f"The equation for the probability is the integral of the probability density:")
    print(f"P = Integral from {x1} to {x2} of 2 * sin^2({n}*pi*x) dx")
    print(f"\nThe calculated probability is: {probability}")
    
    # This result makes physical sense: for n=2, there is a node
    # (zero probability) at the center of the box (x=0.5a).
    # Since the interval [0.495a, 0.505a] is very small and centered
    # on this node, the probability is expected to be very small, but not exactly zero.

if __name__ == "__main__":
    main()
