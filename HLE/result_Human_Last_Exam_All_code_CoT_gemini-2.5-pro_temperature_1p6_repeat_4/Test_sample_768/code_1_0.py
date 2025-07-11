import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function calculates the value of the definite integral.
    """
    # Define the golden ratio
    phi = (1 + np.sqrt(5)) / 2

    # Define the integration limits
    lower_limit = 0
    upper_limit = phi**3 - 1

    # Define the simplified integrand function f(x)
    # f(x) = cos(log(1 + exp(arctan(log(cos(x/e))))))
    def integrand(x):
        # The argument of the innermost log, cos(x/e), is always positive
        # on the integration interval [0, phi^3-1], so the function is well-defined.
        cos_val = np.cos(x / np.e)
        log_cos_val = np.log(cos_val)
        arctan_val = np.arctan(log_cos_val)
        exp_val = np.exp(arctan_val)
        log_val = np.log(1 + exp_val)
        return np.cos(log_val)

    # Perform the numerical integration
    # quad returns the result and an estimated error
    integral_value, integral_error = quad(integrand, lower_limit, upper_limit)

    # As requested, output the numbers in the final equation.
    # The final equation is: Integral from `a` to `b` of f(x)dx = value
    print("The evaluation of the definite integral results in the following equation:")
    print(f"Integral from {lower_limit} to {upper_limit:.10f} = {integral_value:.10f}")
    
solve_integral()