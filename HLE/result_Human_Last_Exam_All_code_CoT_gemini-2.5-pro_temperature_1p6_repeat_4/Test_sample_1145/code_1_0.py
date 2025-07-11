import numpy as np
from scipy import integrate

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box in a specific interval.
    """
    # Define the quantum state and interval boundaries
    n = 2
    x1_factor = 0.495
    x2_factor = 0.505
    
    # We can set the box length 'a' to 1 for simplicity, as it cancels out
    # in the final dimensionless probability.
    a = 1.0
    x1 = x1_factor * a
    x2 = x2_factor * a

    # Define the probability density function, |ψ(x)|²
    # |ψ(x)|² = (2/a) * sin²(nπx / a)
    def probability_density(x, n_val, a_val):
        return (2.0 / a_val) * (np.sin(n_val * np.pi * x / a_val))**2

    # Perform the numerical integration
    # The result is a tuple (integral_value, estimated_error)
    prob, err = integrate.quad(probability_density, x1, x2, args=(n, a))
    
    # Print the equation and the final answer
    print("The probability P of finding the particle in the 1D box is given by the integral of the probability density |ψn(x)|².")
    print(f"For n = {n}, the interval is from {x1_factor}a to {x2_factor}a.")
    print(f"The equation to solve is:")
    print(f"P = ∫ from {x1_factor}a to {x2_factor}a of (2/a) * sin²({n}πx/a) dx")
    print("\nCalculated Probability:")
    print(prob)

solve_probability()
<<<1.315538563388133e-05>>>