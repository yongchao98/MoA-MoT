import numpy as np

def solve():
    """
    This function calculates the analytical approximation for the integral
    I(epsilon) = integral from 0 to 15 of 1 / (epsilon + 9*x^5 + 5*x^6 + 9*x^8) dx
    for small epsilon.
    """
    # The function in the denominator is f(x) = 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0.
    # For small x, the dominant term is a * x^d.
    a = 9.0
    d = 5.0

    # The integral is approximated as I(eps) ~ C * eps^p.
    # The power p is given by p = (1-d)/d.
    p = (1.0 - d) / d

    # The coefficient C is given by C = (1/d) * a**(-1/d) * (pi / sin(pi/d)).
    C = (1.0 / d) * (a**(-1.0 / d)) * (np.pi / np.sin(np.pi / d))

    # Print the final derived formula with the calculated numerical values.
    # The problem asks to output each number in the final equation.
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(epsilon) \u2248 {C:.4f} * epsilon^({p})")
    
    # We can also express C and p as fractions for clarity
    print("\nIn fractional form, the formula is:")
    print(f"I(epsilon) \u2248 C * epsilon^(-4/5)")
    print(f"where the coefficient C = (\u03C0 / (5 * sin(\u03C0/5))) * 9^(-1/5) \u2248 {C:.4f}")

solve()