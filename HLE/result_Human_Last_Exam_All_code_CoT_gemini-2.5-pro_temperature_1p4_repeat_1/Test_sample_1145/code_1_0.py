import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a given interval
    for a specified quantum state in a 1D box.
    """
    # Define the constants for the problem
    n = 2
    # We can set the length of the box 'a' to 1.0 for simplicity in calculation,
    # as the probability is a dimensionless quantity.
    a = 1.0
    x1 = 0.495 * a
    x2 = 0.505 * a

    # The integrand is the probability density function, |psi(x)|^2.
    # |psi_n(x)|^2 = (2/a) * sin^2(n*pi*x/a)
    def probability_density(x, n_val, a_val):
        return (2 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # Perform the numerical integration using scipy.integrate.quad
    # The quad function returns the integral result and an estimated error.
    probability, _ = quad(probability_density, x1, x2, args=(n, a))

    # Print the explanation and the final equation with all numbers
    print("The probability P is found by integrating the probability density function for n=2")
    print("from x1 = 0.495a to x2 = 0.505a.")
    print("\nThe equation being solved is:")
    print(f"P = \u222B[{x1:.3f}a, {x2:.3f}a] (2/a) * sin\u00b2({n}\u03c0x/a) dx")
    print(f"\nCalculated Probability: {probability}")
    
    # Return the final answer in the specified format
    return probability

# Run the function and capture the result for the final answer format
final_probability = solve_probability()
print(f"<<<{final_probability}>>>")
