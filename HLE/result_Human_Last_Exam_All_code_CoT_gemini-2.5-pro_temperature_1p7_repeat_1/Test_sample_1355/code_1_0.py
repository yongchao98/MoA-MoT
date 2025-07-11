import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

def calculate_moment_ratio():
    """
    This function calculates the ratio of the fourth statistical moment of the dimensionless
    conductance to its average value for a critical disordered Majorana wire.

    The conductance 'g' is given by g = 1 / (1 + x), where x follows a Gamma distribution
    P(x) = (1/sqrt(2*pi*x)) * exp(-x/2).

    The n-th moment <g^n> is given by the integral:
    <g^n> = Integral from 0 to infinity of (1+x)^-n * P(x) dx.
    """

    # Probability density function for the random variable lambda
    def p(x):
        if x <= 0:
            return 0
        return (1 / np.sqrt(2 * np.pi * x)) * np.exp(-x / 2)

    # Integrand for the n-th moment of g
    def integrand(x, n):
        return (1 + x)**(-n) * p(x)

    # Calculate the average conductance <g> (n=1)
    g_avg, _ = quad(integrand, 0, np.inf, args=(1,))

    # Calculate the fourth moment <g^4> (n=4)
    g4_moment, _ = quad(integrand, 0, np.inf, args=(4,))
    
    # Calculate the ratio
    ratio = g4_moment / g_avg

    print("Equation: ratio = <g^4> / <g>")
    print(f"The fourth moment of the conductance is: <g^4> = {g4_moment}")
    print(f"The average value of the conductance is: <g> = {g_avg}")
    print(f"The ratio is: {ratio}")

calculate_moment_ratio()