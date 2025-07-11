import numpy as np
from scipy.integrate import quad

def evaluate_integral():
    """
    Evaluates the definite integral specified by the user.
    The integral is I = integral from 0 to phi^3-1 of Re[1 / (1 + exp(arctan(ln(cos(x/e)))))^i] dx
    """

    # Step 1: Define the constants and integration limits
    phi = (1 + np.sqrt(5)) / 2
    e_val = np.e
    
    lower_bound = 0
    upper_bound = phi**3 - 1

    # Step 2: Define the integrand function
    # As derived, the integrand simplifies to: cos(ln(1 + exp(arctan(ln(cos(x/e))))))
    def integrand(x):
        # To avoid potential domain errors with edge cases, though not expected in this interval
        cos_x_e = np.cos(x / e_val)
        if cos_x_e <= 0:
            return 0
        
        # Calculate the simplified expression
        val = np.cos(np.log(1 + np.exp(np.arctan(np.log(cos_x_e)))))
        return val

    # Step 3: Perform the numerical integration
    result, error = quad(integrand, lower_bound, upper_bound)

    # Step 4: Print the numbers in the final equation
    print(f"The final equation is: integral from {lower_bound} to {upper_bound} of f(x)dx = {result}")
    print("\nWhere the numbers in the equation are:")
    print(f"Lower bound: {lower_bound}")
    print(f"Upper bound (phi^3 - 1): {upper_bound}")
    print(f"Value of the integral: {result}")
    
evaluate_integral()