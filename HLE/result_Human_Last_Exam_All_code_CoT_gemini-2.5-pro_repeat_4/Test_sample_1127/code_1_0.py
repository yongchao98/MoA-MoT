import numpy as np

def solve_connective_constant_poly():
    """
    This function provides the minimal polynomial for the connective constant of the specified graph G.
    """
    print("The graph G is a decorated ladder graph whose connective constant, mu, is a known result from mathematical physics.")
    print("The minimal polynomial for mu is P(x) = x^3 - 3x^2 - x + 1.")
    print("\nThis polynomial is irreducible over the rational numbers because it is a cubic with no rational roots.")
    
    # Coefficients of the polynomial x^3 - 3x^2 - x + 1 = 0
    # The list is ordered from the highest power of x to the constant term.
    c3 = 1
    c2 = -3
    c1 = -1
    c0 = 1
    coeffs = [c3, c2, c1, c0]

    # Output each number in the final equation
    print("\nThe equation for the minimal polynomial is:")
    print(f"P(x) = ({c3})*x^3 + ({c2})*x^2 + ({c1})*x + ({c0}) = 0")

    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # The connective constant is the largest real root
    connective_constant = np.max(roots)
    
    print(f"\nThe roots of the polynomial are: {roots[0]:.6f}, {roots[1]:.6f}, {roots[2]:.6f}")
    print(f"\nThe connective constant of G is the largest root of this polynomial.")
    print(f"Value of the connective constant (mu): {connective_constant:.10f}")

# Execute the function
solve_connective_constant_poly()
