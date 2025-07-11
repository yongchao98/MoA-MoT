import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function evaluates the definite integral:
    ∫[0 to φ³-1] Re[1 / (1 + exp(arctan(ln(cos(x/e)))))**i] dx
    """

    # Step 1: Define constants and integration limits
    phi = (1 + np.sqrt(5)) / 2
    e = np.e
    
    lower_limit = 0
    # The upper limit is φ³-1.
    # φ² = φ + 1
    # φ³ = φ² * φ = (φ+1) * φ = φ² + φ = (φ+1) + φ = 2φ + 1
    # So, φ³ - 1 = 2φ = 1 + sqrt(5)
    upper_limit = 2 * phi

    # Step 2: Define the integrand function
    # The integrand simplifies to cos(ln(1 + exp(arctan(ln(cos(x/e))))))
    def integrand(x):
        # The argument x/e is always within [-π/2, π/2] for the integration range,
        # so cos(x/e) is positive and the expression is well-defined.
        cos_val = np.cos(x / e)
        log_cos_val = np.log(cos_val)
        arctan_val = np.arctan(log_cos_val)
        exp_val = np.exp(arctan_val)
        log_val = np.log(1 + exp_val)
        return np.cos(log_val)

    # Step 3: Perform numerical integration
    result, error = quad(integrand, lower_limit, upper_limit)

    # Step 4: Output the results as a final equation
    print("The integral evaluation involves the following numbers:")
    print(f"Lower Limit = {lower_limit}")
    print(f"Upper Limit (φ³ - 1) = {upper_limit}")
    print("\nThe final equation is:")
    print(f"Integral from {lower_limit} to {upper_limit:.5f} = {result:.5f}")

solve_integral()