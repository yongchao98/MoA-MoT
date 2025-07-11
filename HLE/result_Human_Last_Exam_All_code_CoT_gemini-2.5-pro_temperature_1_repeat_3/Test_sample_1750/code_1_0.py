import numpy as np
from scipy.integrate import quad
from fractions import Fraction

def solve_integral():
    """
    This function solves the definite integral numerically and returns the
    answer as a fraction of pi.
    """
    
    # The integral simplifies to I = integral from 0 to pi of (sin(2x)(2cos(2x)-1))^50 dx
    # Let's perform a substitution u = 2x, so dx = du/2.
    # The limits for x from 0 to pi become u from 0 to 2*pi.
    # I = (1/2) * integral from 0 to 2*pi of (sin(u)(2cos(u)-1))^50 du.
    # The integrand f(u) = (sin(u)(2cos(u)-1))^50 has a period of 2*pi.
    # Let's check the symmetry f(2*pi - u) vs f(u).
    # f(2*pi - u) = (sin(2*pi-u)(2cos(2*pi-u)-1))^50 = (-sin(u)(2cos(u)-1))^50 = f(u).
    # This symmetry means integral from 0 to 2*pi is 2 * integral from 0 to pi.
    # So, I = (1/2) * 2 * integral from 0 to pi of (sin(u)(2cos(u)-1))^50 du
    # I = integral from 0 to pi of (sin(u)(2cos(u)-1))^50 du.
    
    def integrand(u):
        return (np.sin(u) * (2 * np.cos(u) - 1))**50

    # Perform numerical integration
    integral_value, error = quad(integrand, 0, np.pi)

    # The result is expected to be a rational multiple of pi.
    # Let's find the ratio.
    ratio_to_pi = integral_value / np.pi
    
    # Convert the ratio to a simple fraction.
    # The limit_denominator is set to a large value to handle potential
    # floating point inaccuracies, but for such problems, the true
    # denominator is usually small.
    fraction = Fraction(ratio_to_pi).limit_denominator(1000000)
    
    numerator = fraction.numerator
    denominator = fraction.denominator

    print(f"The simplified integral is integral from 0 to pi of (sin(x)*(2*cos(x)-1))^50 dx.")
    print(f"Numerical integration result: {integral_value}")
    print(f"Result divided by pi: {ratio_to_pi}")
    print(f"The fraction of pi is: {numerator} / {denominator}")
    print(f"The final answer is: pi * {numerator} / {denominator}")

solve_integral()
