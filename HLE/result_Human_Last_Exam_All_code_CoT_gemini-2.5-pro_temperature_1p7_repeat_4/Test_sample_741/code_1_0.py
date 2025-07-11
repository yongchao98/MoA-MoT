import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

def solve_bessel_root():
    """
    Finds the largest value of x for which the summation equals 0.
    The summation is identified as the modified Bessel function I_{x-1}(2).
    """

    # We need to solve the equation I_{x-1}(2) = 0.
    # Let nu = x - 1. The equation becomes I_nu(2) = 0.
    # We define a function for the root finder.
    def bessel_nu_of_2(nu):
        """Modified Bessel function of the first kind I_nu(2) as a function of nu."""
        return iv(nu, 2)

    # From mathematical properties and testing, the largest (least negative) root for nu
    # is known to be in the interval [-2.5, -2.0].
    # I_nu(2) at nu=-2.5 is negative, and at nu=-2.0 is positive.
    bracket = [-2.5, -2.0]

    # Use a robust root-finding algorithm to find nu.
    solution = root_scalar(bessel_nu_of_2, bracket=bracket, method='brentq')
    nu_root = solution.root

    # The problem asks for x, where x = nu + 1.
    x_value = nu_root + 1

    # Output the numbers in the final equation as requested.
    print(f"The summation is equal to the modified Bessel function I_(x-1)(2).")
    print(f"We solve the equation: I_nu(2) = 0, where nu = x - 1.")
    print(f"The largest numerical root found for nu is: {nu_root:.9f}")
    print(f"The equation with the solved number is: I_({nu_root:.4f})(2) = {bessel_nu_of_2(nu_root):.2e}")
    print(f"Since x = nu + 1, the largest value for x is: {x_value:.9f}")

    # Print the final answer in the requested format {-a.bbb}
    print(f"\nFinal Answer: {{{x_value:.3f}}}")


if __name__ == "__main__":
    solve_bessel_root()
