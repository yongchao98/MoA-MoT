import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Solves the specified definite integral numerically.
    """
    # The golden ratio, phi
    phi = (1 + np.sqrt(5)) / 2

    # The upper limit of integration is phi^3 - 1
    # We can simplify this using the property phi^2 = phi + 1:
    # phi^3 = phi * phi^2 = phi * (phi + 1) = phi^2 + phi = (phi + 1) + phi = 2*phi + 1
    # So, phi^3 - 1 = 2*phi
    upper_limit = 2 * phi

    # The numbers from the final equation for the integral can be considered the limits
    # and the result itself.
    lower_limit = 0

    print(f"Evaluating the integral from {lower_limit} to {upper_limit:.10f}")

    def integrand(x):
        """
        This is the simplified real-valued function to be integrated.
        integrand(x) = cos(ln(1 + exp(arctan(ln(cos(x/e))))))
        """
        # To avoid math domain errors for values outside the integration range,
        # we check if cos(x/e) is positive. It is for our integration interval.
        cos_val = np.cos(x / np.e)
        
        log_cos_val = np.log(cos_val)
        arctan_val = np.arctan(log_cos_val)
        exp_val = np.exp(arctan_val)
        
        base_for_log = 1 + exp_val
        final_arg_for_cos = np.log(base_for_log)
        
        return np.cos(final_arg_for_cos)

    # Perform the numerical integration using scipy.integrate.quad
    result, error = quad(integrand, lower_limit, upper_limit)

    print(f"The equation is: Integral from {lower_limit} to (phi^3 - 1) is equal to the result.")
    print(f"Numbers from the equation are:")
    print(f"Lower Limit: {lower_limit}")
    print(f"Upper Limit (phi^3 - 1): {upper_limit}")
    print(f"Result: {result}")
    
    # Returning the final numerical result for the 'answer' block.
    return result

# Run the solver and print the final answer in the required format
final_result = solve_integral()
print(f"Final Answer: {final_result}")
