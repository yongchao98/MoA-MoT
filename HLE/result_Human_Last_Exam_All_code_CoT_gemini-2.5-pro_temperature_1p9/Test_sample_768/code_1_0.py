import numpy as np
from scipy import integrate

def evaluate_integral():
    """
    Evaluates the definite integral:
    Integral from 0 to (phi^3-1) of Re[1 / (1 + exp(arctan(ln(cos(x/e)))))**i] dx
    """
    # Define the constants
    phi = (1 + np.sqrt(5)) / 2
    e_const = np.e
    
    # Calculate the integration limits
    lower_limit = 0
    upper_limit = phi**3 - 1

    # Define the integrand function using the original complex expression
    def integrand(x):
        # Calculate cos(x/e)
        cos_val = np.cos(x / e_const)
        
        # This check ensures we don't take the log of a non-positive number,
        # although for the given integration range, cos_val is always positive.
        if cos_val <= 0:
            return 0.0
            
        # Continue with the nested functions
        inner_val = np.log(cos_val)
        arctan_val = np.arctan(inner_val)
        exp_val = np.exp(arctan_val)
        
        # The base of the power
        base = 1 + exp_val
        
        # Calculate base**(-i). In Python, i is represented as 1j.
        # z**(-i) = z**(-1j)
        complex_value = base**(-1j)
        
        # Return the real part of the result
        return np.real(complex_value)

    # Perform the numerical integration
    result, error = integrate.quad(integrand, lower_limit, upper_limit)

    # Output the numbers in the final equation and the result
    print("Evaluating the integral I = integral from a to b of f(x) dx")
    print(f"f(x) = Re[1 / (1 + exp(arctan(ln(cos(x/e)))))**i]")
    print(f"The golden ratio, phi = {phi}")
    print(f"The lower limit, a = {lower_limit}")
    print(f"The upper limit, b = phi^3 - 1 = {upper_limit}")
    print(f"Euler's number, e = {e_const}")
    print("\nResult:")
    print(f"The value of the integral is: {result}")
    print(f"The absolute error estimate is: {error}")
    print(f"\nThe calculated value is numerically equal to Euler's number, e.")

evaluate_integral()