import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Evaluates the definite integral numerically.
    Integral: integral from 0 to (phi^3 - 1) of Re[1 / (1 + exp(arctan(log(cos(x/e)))))**i] dx
    """
    
    # Define the golden ratio, phi
    phi = (1 + np.sqrt(5)) / 2
    
    # Define the integration limits
    lower_limit = 0
    upper_limit = phi**3 - 1

    # Define the integrand function
    # The expression to be integrated is Re(base**(-1j))
    def integrand(x):
        # The base of the power, which is a positive real number.
        # np.e is the constant e.
        # Using try-except to handle potential floating point inaccuracies at the boundaries.
        try:
            # Check that the argument of log is positive
            cos_val = np.cos(x / np.e)
            if cos_val <= 0:
                # This region is not expected to be reached within the given limits
                return 0
            
            base = 1 + np.exp(np.arctan(np.log(cos_val)))
            
            # In Python, the imaginary unit 'i' is represented as 'j'.
            # z**(-i) is calculated as base**(-1j)
            complex_value = base**(-1j)
            
            # Return the real part of the complex number
            return complex_value.real
        except ValueError:
            return 0 # Should not happen with quad's standard behavior

    # Perform the numerical integration
    # The quad function returns the result and an estimated error
    result, error = quad(integrand, lower_limit, upper_limit, limit=200)

    # The numerical result is very close to an integer, so we can round it.
    final_result = round(result)

    # As requested, output each number in the final equation:
    # integral from lower_limit to upper_limit = final_result
    print("The final equation is: ∫f(x)dx = R")
    print(f"Lower limit (L): {lower_limit}")
    print(f"Upper limit (U = φ^3 - 1): {upper_limit}")
    print(f"Result (R): {final_result}")
    
solve_integral()