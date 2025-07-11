import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Evaluates the definite integral using numerical methods.
    """
    
    # Define the golden ratio
    phi = (1 + np.sqrt(5)) / 2
    
    # Define the integration limits
    a = 0
    b = phi**3 - 1
    
    # Define the simplified real-valued integrand
    # f(x) = cos(ln(1 + exp(atan(ln(cos(x/e))))))
    def integrand(x):
        # Handle the edge case x=0 where cos(0)=1, ln(1)=0
        if x == 0:
            return np.cos(np.log(2))
        
        inner_val = np.cos(x / np.e)
        # Avoid math domain error for inner_val <= 0, though not expected in this range
        if inner_val <= 0:
            return 0

        arg_of_cos = np.log(1 + np.exp(np.arctan(np.log(inner_val))))
        return np.cos(arg_of_cos)

    # Perform the numerical integration
    result, error = quad(integrand, a, b)
    
    # The user asked me to be a helpful assistant, so I'm not only solving the integral but also presenting the inputs.
    print(f"Golden ratio (φ): {phi}")
    print(f"Upper limit (φ^3 - 1): {b}")
    print(f"Value of the integral: {result}")

solve_integral()
