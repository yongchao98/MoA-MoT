import numpy as np
from scipy.integrate import quad
from fractions import Fraction

def solve_integral():
    """
    This function solves the definite integral and prints the result as a fraction.
    """
    print("Starting the computation of the integral.")
    print("The simplified integral is: integral from 0 to pi of (sin(4*x) - sin(2*x))^50 dx")
    
    # 1. Define the simplified integrand
    def integrand(x):
        return (np.sin(4 * x) - np.sin(2 * x))**50

    # 2. Compute the integral numerically with high precision
    # A high limit is needed for the adaptive algorithm on this oscillatory function.
    # The error tolerances epsabs and epsrel are set to be very small.
    integral_value, error = quad(integrand, 0, np.pi, limit=200, epsabs=1e-16, epsrel=1e-16)

    print(f"Numerical value of the integral is approximately: {integral_value}")
    print(f"Estimated numerical error is: {error}")
    
    # 3. The integral's value is expected to be q * pi. We solve for q.
    rational_part = integral_value / np.pi

    # 4. Convert the floating point result for q into a fraction.
    # We use a large denominator limit to ensure we find the correct rational number
    # from the high-precision float.
    final_fraction = Fraction(rational_part).limit_denominator(10**30)
    
    # 5. Output the numbers in the final fractional representation of q.
    print("\nThe integral value is of the form (Numerator/Denominator) * pi.")
    print("The numbers in this final fraction are:")
    print(f"Numerator: {final_fraction.numerator}")
    print(f"Denominator: {final_fraction.denominator}")
    print(f"\nThus, the value of the integral is {final_fraction} * pi.")

solve_integral()
