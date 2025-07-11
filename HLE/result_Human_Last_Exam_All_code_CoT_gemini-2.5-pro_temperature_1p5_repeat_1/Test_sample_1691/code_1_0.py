import numpy as np

def solve_integral_approximation():
    """
    This function develops and prints an analytical formula that approximates
    the integral I(epsilon) for small epsilon.
    """
    
    # The integral is I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9x^5 + 5x^6 + 9x^8) dx.
    # For small x, the denominator is dominated by epsilon + a*x^n.
    a = 9.0
    n = 5.0
    
    # The general approximation for small epsilon is I(epsilon) ~ C * epsilon^p, where
    # p = -(n-1)/n
    # C = pi / (n * a^(1/n) * sin(pi/n))
    
    # Calculate the power p
    p = -(n - 1.0) / n
    
    # Calculate the coefficient C
    C = np.pi / (n * np.power(a, 1.0/n) * np.sin(np.pi/n))
    
    print("For small epsilon, the integral I(epsilon) can be approximated by the formula:")
    print("I(epsilon) ≈ C * ε^p\n")
    print("This approximation is derived by considering the dominant term in the denominator near x=0, which is 9.0 * x^5.0.\n")
    
    print("The power 'p' is calculated as:")
    print(f"p = -(n - 1) / n = -({n:.1f} - 1.0) / {n:.1f} = {p:.1f}")
    
    print("\nThe coefficient 'C' is calculated using the formula:")
    print(f"C = π / (n * a^(1/n) * sin(π/n))")
    print("Plugging in the values n=5.0 and a=9.0:")
    print(f"C = π / ({n:.1f} * {a:.1f}^(1/{n:.1f}) * sin(π/{n:.1f}))")
    print(f"C ≈ {np.pi:.6f} / ({n:.1f} * {np.power(a, 1.0/n):.6f} * {np.sin(np.pi/n):.6f})")
    print(f"C ≈ {C:.6f}\n")
    
    print("Therefore, the final analytical approximation is:")
    # Using 'ε' for epsilon for better readability
    print(f"I(ε) ≈ {C:.6f} * ε^({p:.1f})")

solve_integral_approximation()