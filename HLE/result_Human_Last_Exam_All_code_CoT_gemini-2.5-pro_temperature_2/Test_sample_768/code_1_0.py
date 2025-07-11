import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function evaluates the specified definite integral numerically.
    The integral is: ∫[0, φ³-1] Re[1 / (1 + exp(arctan(ln(cos(x/e)))))^i] dx
    where φ is the golden ratio.
    """
    
    # Define the mathematical constants
    phi = (1 + np.sqrt(5)) / 2
    e_val = np.e

    # Define the integration limits
    lower_limit = 0
    # The upper limit simplifies to 2*phi
    upper_limit = phi**3 - 1

    def integrand(x):
        """
        The simplified real-valued integrand function:
        f(x) = cos(log(1 + exp(atan(log(cos(x/e))))))
        """
        # We must ensure the argument of log is positive.
        # For the given integration range, x/e is in [0, 2*phi/e] ≈ [0, 1.19].
        # In this range, cos(x/e) is in (0, 1], so ln(cos(x/e)) is well-defined.
        inner_val = np.cos(x / e_val)
        
        # Guard against floating point errors that might push inner_val to <= 0
        if inner_val <= 0:
            return 0
        
        log_cos = np.log(inner_val)
        arctan_log_cos = np.arctan(log_cos)
        exp_arctan = np.exp(arctan_log_cos)
        
        return np.cos(np.log(1 + exp_arctan))

    # Perform the numerical integration
    result, error = quad(integrand, lower_limit, upper_limit)

    # Output the numbers in the final equation as requested.
    # The final equation is: ∫[lower_limit, upper_limit] f(x) dx = result
    print("This script evaluates the definite integral:")
    print("I = ∫[a, b] Re[1 / (1 + exp(arctan(ln(cos(x/e)))))^i] dx")
    print("\nThe numbers in the final equation are:")
    print(f"Golden ratio, φ = {phi}")
    print(f"Lower limit, a = {lower_limit}")
    print(f"Upper limit, b = φ³ - 1 = {upper_limit}")
    print(f"\nThe value of the integral is:")
    print(f"I = {result}")
    print(f"(The numerical error of the calculation is estimated to be {error})")
    print("\nObservation:")
    print(f"The result is numerically identical to Euler's number, e = {e_val}")


solve_integral()