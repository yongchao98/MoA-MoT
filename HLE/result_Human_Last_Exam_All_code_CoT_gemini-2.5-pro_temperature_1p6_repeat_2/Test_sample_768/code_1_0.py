import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function evaluates the definite integral
    ∫[0, φ^3-1] Re[1 / (1 + exp(arctan(ln(cos(x/e)))))**i] dx
    by first simplifying the integrand and then using numerical integration.
    """
    # 1. Simplify the integrand as derived in the thinking process.
    def integrand(x):
        """
        The simplified integrand: cos(ln(1 + exp(arctan(ln(cos(x/e)))))).
        """
        # Using np.e for Euler's number
        val_e = np.e
        
        # Check domain for ln, although for the given limits, cos(...) > 0.
        cos_val = np.cos(x / val_e)
        if cos_val <= 0:
            return 0 # Should not be reached
        
        # Sequentially calculate the parts of the expression
        log_of_cos = np.log(cos_val)
        arctan_of_log = np.arctan(log_of_cos)
        exp_of_arctan = np.exp(arctan_of_log)
        
        # This is the 'Z' from the derivation
        Z = 1 + exp_of_arctan
        
        # The final result of the integrand is cos(ln(Z))
        return np.cos(np.log(Z))

    # 2. Define the constants and limits of integration.
    # Golden Ratio, φ
    phi = (1 + np.sqrt(5)) / 2
    
    # Upper limit of integration is φ^3 - 1
    # Note: φ^3 - 1 = (2*φ + 1) - 1 = 2*φ = 1 + sqrt(5)
    upper_limit = phi**3 - 1
    lower_limit = 0

    # 3. Perform numerical integration.
    # The quad function returns the integral result and an error estimate.
    integral_value, _ = quad(integrand, lower_limit, upper_limit)

    # 4. Output the result in the requested format.
    # "Final equation" will be Integral = [value]
    # "each number" will be the constants and the final value.
    
    print(f"The value of the golden ratio is φ = {phi}")
    print(f"The upper limit of the integral is φ^3 - 1 = {upper_limit}")
    
    # Final equation and its single number (the result)
    print("\nThe final equation is:")
    print(f"I = {integral_value}")

# Execute the function to solve and print the result.
solve_integral()