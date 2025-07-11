import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Evaluates the definite integral numerically.
    """
    # Define the golden ratio, phi
    phi = (1 + np.sqrt(5)) / 2

    # Define the upper limit of integration
    # Note: phi**3 - 1 = (2*phi + 1) - 1 = 2*phi
    lower_limit = 0
    upper_limit = phi**3 - 1

    def integrand(x):
        """
        The simplified real-valued integrand function.
        f(x) = cos(ln(1 + exp(arctan(ln(cos(x/e))))))
        """
        # Python's built-in math constants are used for precision.
        # Note: np.e is an acceptable alias for math.e
        inner_val = np.cos(x / np.e)
        
        # Add a small epsilon to prevent log(0) for arguments at the edge of validity,
        # although for the given integration range, cos(x/e) is strictly positive.
        epsilon = 1e-15
        log_val = np.log(inner_val + epsilon)
        
        arctan_val = np.arctan(log_val)
        exp_val = np.exp(arctan_val)
        
        final_log = np.log(1 + exp_val)
        
        return np.cos(final_log)

    # Perform the numerical integration using scipy.integrate.quad
    result, error = quad(integrand, lower_limit, upper_limit, limit=200)

    # The problem asks to output the numbers in the final equation.
    # The final equation is: integral from 0 to (phi^3-1) = result.
    print(f"The integral is evaluated from {lower_limit} to {upper_limit:.10f}")
    print(f"The calculated value of the integral is: {result:.10f}")
    
    # The numerical result is remarkably close to phi^2.
    # We can print phi^2 for comparison.
    # print(f"For comparison, the value of phi^2 is: {phi**2:.10f}")
    
solve_integral()