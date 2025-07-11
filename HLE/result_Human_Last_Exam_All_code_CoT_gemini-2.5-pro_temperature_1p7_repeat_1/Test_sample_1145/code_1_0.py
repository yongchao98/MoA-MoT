import numpy as np
from scipy.integrate import quad

def solve_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    in a given interval for a specific quantum state n.
    """
    # Define the quantum state
    n = 2

    # Define the integration limits as a fraction of the box length 'a'
    # u1 = x1/a and u2 = x2/a
    u1 = 0.495
    u2 = 0.505

    # Define the probability density function for a box of length a=1
    # P(u) = 2 * sin^2(n * pi * u)
    def probability_density(u, n_level):
        return 2 * (np.sin(n_level * np.pi * u))**2

    # Perform the numerical integration using scipy.integrate.quad
    # The quad function returns the result of the integral and an error estimate.
    probability, _ = quad(probability_density, u1, u2, args=(n,))

    # Print the equation with the specific values used for the calculation
    print("To find the probability, we solve the following integral:")
    print(f"P = Integral from u1={u1} to u2={u2} of [2 * sin^2({n} * pi * u)] du")
    
    # Print the final calculated probability
    print("\nThe result of the integration is:")
    print(f"Probability = {probability}")

solve_probability()
<<<6.579736113481238e-06>>>