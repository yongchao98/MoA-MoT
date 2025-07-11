import numpy as np

def solve_integral_approximation():
    """
    This function derives and prints the analytical approximation for the integral
    I(epsilon) for small epsilon.
    """
    
    # The integral is I(epsilon) = integral_0^15 dx / (epsilon + 9x^5 + 5x^6 + 9x^8)
    # For small x, the denominator is dominated by the term with the lowest power of x.
    # Denominator approx epsilon + c0 * x^n
    c0 = 9.0
    n = 5.0
    
    # The approximation for the integral is I(epsilon) ~ C * epsilon^alpha
    # The exponent alpha is given by -(n-1)/n
    alpha_num = -(n - 1.0)
    alpha_den = n
    alpha = alpha_num / alpha_den
    
    # The coefficient C is given by (c0^(-1/n) * pi) / (n * sin(pi/n))
    # This comes from evaluating the integral: integral_0^inf dx / (epsilon + c0*x^n)
    C = (c0**(-1.0/n) * np.pi) / (n * np.sin(np.pi/n))
    
    print("The analytical formula that approximates I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) = C * epsilon^alpha")
    
    print("\nBased on the leading term 9.0 * x^5.0 in the denominator, the constants are derived as follows:")
    print(f"alpha = -({n} - 1.0) / {n}")
    print(f"C = ({c0}^(-1.0/{n}) * pi) / ({n} * sin(pi/{n}))")
    
    print("\nSubstituting the numerical values into these expressions:")
    print(f"alpha = {alpha}")
    print(f"C approx {C:.6f}")
    
    print("\nThus, the final approximate analytical formula is:")
    print(f"I(epsilon) approx ({C:.6f}) * epsilon^({alpha})")
    
    print("\nThe full analytical expression with each number is:")
    print(f"I(epsilon) = ({c0}^(-1.0/{n}) * pi / ({n} * sin(pi/{n}))) * epsilon^(-({n}-1.0)/{n})")

solve_integral_approximation()