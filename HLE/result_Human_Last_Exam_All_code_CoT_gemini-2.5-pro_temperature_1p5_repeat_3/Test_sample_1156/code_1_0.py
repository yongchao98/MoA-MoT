import numpy as np
from scipy.integrate import quad

def main():
    """
    Calculates and prints the parameters for the normalised invariant density
    for the map T(x) = 1/sqrt(x) mod 1.
    """
    
    # The un-normalised invariant density function is h(x) = 1/sqrt(x) - 1/(1+x).
    # We need to handle the singularity at x=0 for the numerical integration.
    # The integral is well-defined.
    def h(x):
        if x == 0:
            return np.inf
        return 1 / np.sqrt(x) - 1 / (1 + x)

    # Calculate the integral of h(x) from 0 to 1 to find the normalization constant.
    # The integral of 1/sqrt(x) is 2*sqrt(x).
    # The integral of 1/(1+x) is ln(1+x).
    # So, integral from 0 to 1 is [2*sqrt(x) - ln(1+x)] from 0 to 1
    # = (2*sqrt(1) - ln(2)) - (2*sqrt(0) - ln(1)) = 2 - ln(2).
    # We use scipy.integrate.quad to verify this numerically.
    integral_value, error = quad(h, 0, 1)

    # The normalised density is rho(x) = (1/integral_value) * h(x).
    # The equation for the density is of the form:
    # rho(x) = (1 / (p1 - ln(p2))) * (1/x**p3 - 1/(p4 + x))
    
    p1 = 2.0
    p2 = 2.0
    p3 = 0.5
    p4 = 1.0

    print("The normalised density of the invariant measure is rho(x), where:")
    print("rho(x) = C * (1/x**p3 - 1/(p4 + x))")
    print("The parameters for this equation are:")
    print(f"p1 (in constant C) = {p1}")
    print(f"p2 (in constant C) = {p2}")
    print(f"p3 (exponent) = {p3}")
    print(f"p4 (in denominator) = {p4}")
    
    print("\nThe normalization constant C is given by C = 1 / (p1 - ln(p2))")
    C = 1 / integral_value
    print(f"Numerically, C = 1 / {integral_value:.8f} = {C:.8f}")
    
    analytical_C = 1 / (p1 - np.log(p2))
    print(f"Analytically, C = 1 / (2 - ln(2)) = {analytical_C:.8f}")

if __name__ == '__main__':
    main()
