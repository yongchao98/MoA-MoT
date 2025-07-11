import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Solves the definite integral by first simplifying the integrand and then performing numerical integration.
    """
    # Step 1 & 2: Define constants and integration limits
    golden_ratio = (1 + np.sqrt(5)) / 2
    lower_limit = 0
    upper_limit = golden_ratio**3 - 1

    # Step 3: Define the integrand function based on the simplification
    def integrand(x):
        """
        This is the simplified real part of the original complex expression.
        f(x) = cos(log(1 + exp(arctan(log(cos(x/e))))))
        """
        # We checked that cos(x/np.e) is positive over the integration interval,
        # so all intermediate operations are well-defined with real numbers.
        cos_val = np.cos(x / np.e)
        log_cos_val = np.log(cos_val)
        arctan_val = np.arctan(log_cos_val)
        exp_val = np.exp(arctan_val)
        
        # This is the base of the power B(x) from the explanation
        base = 1 + exp_val
        
        # The real part is cos(ln(B(x)))
        return np.cos(np.log(base))

    # Step 4: Perform the numerical integration
    integral_value, error = quad(integrand, lower_limit, upper_limit)

    # Step 5: Output the numbers in the final equation as requested
    print("The final equation is of the form: Integral = Value")
    print(f"The value of the golden ratio phi is: {golden_ratio}")
    print(f"The upper limit of integration (phi^3 - 1) is: {upper_limit}")
    print(f"The calculated value of the integral is: {integral_value}")

solve_integral()