import numpy as np
from scipy.integrate import quad

def calculate_integral():
    """
    Calculates the value of the definite integral:
    ∫[0, 2] (2**(-1/16) * tan(asin(x**4/(16*sqrt(2)))) + 2**(1/16) * (sin(atan(x/2)))**(1/4)) dx
    """

    # The integrand can be simplified for easier computation.
    # Term 1: 2**(-1/16) * tan(asin(x**4 / (16*sqrt(2))))
    # Let theta = asin(x**4 / (16*sqrt(2))).
    # tan(theta) simplifies to x**4 / sqrt(512 - x**8).
    # So, Term 1 = 2**(-1/16) * x**4 / sqrt(512 - x**8)
    
    # Term 2: 2**(1/16) * (sin(atan(x/2)))**(1/4)
    # Let phi = atan(x/2).
    # sin(phi) simplifies to x / sqrt(x**2 + 4).
    # So, Term 2 = 2**(1/16) * (x / sqrt(x**2 + 4))**(1/4)
    # which is 2**(1/16) * x**(1/4) / (x**2 + 4)**(1/8)

    def integrand(x):
        """The simplified integrand function."""
        # Note: x=0 is a special case. x**(1/4) is 0, and x**4 is 0.
        # The quad function handles endpoints correctly, but we must be careful if writing the code differently.
        if x == 0:
            return 0
            
        term1 = 2**(-1/16) * (x**4) / np.sqrt(512 - x**8)
        term2 = 2**(1/16) * (x**(1/4)) / np.power(x**2 + 4, 1/8)
        
        return term1 + term2

    # Integrate the function from 0 to 2
    # The first number is the constant from the first term in the original expression.
    c1 = 2**(-1/16)
    # The second number is the constant from the second term.
    c2 = 2**(1/16)
    
    integral_value, error_estimate = quad(integrand, 0, 2)
    
    # Per the instruction to output each number in the final equation,
    # we present the final result as a clear statement.
    # Final Equation: ∫[0, 2] (c1 * f1(x) + c2 * f2(x)) dx = integral_value
    print(f"The first constant is c1 = 2**(-1/16) = {c1}")
    print(f"The second constant is c2 = 2**(1/16) = {c2}")
    print(f"The value of the definite integral is: {integral_value}")
    print(f"Estimated error of the calculation is: {error_estimate}")

if __name__ == '__main__':
    calculate_integral()
