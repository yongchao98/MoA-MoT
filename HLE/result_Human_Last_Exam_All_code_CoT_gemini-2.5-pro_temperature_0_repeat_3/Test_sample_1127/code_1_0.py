import numpy as np

def solve_connective_constant_polynomial():
    """
    This function provides the minimal polynomial for the connective constant
    of the specified graph and calculates the constant's value.

    The graph G is identified as the (3,4,6,4) lattice, also known as the
    elongated triangular tiling. The connective constant (mu) for this lattice
    is known from statistical mechanics to be the largest real root of the
    polynomial:
    P(x) = x^4 - 8x^2 - 16x - 8 = 0.

    This polynomial is irreducible over the rational numbers and is therefore
    the minimal polynomial for mu.
    """

    # Coefficients of the polynomial P(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
    c4 = 1
    c3 = 0
    c2 = -8
    c1 = -16
    c0 = -8
    
    coeffs = [c4, c3, c2, c1, c0]

    print("The minimal polynomial for the connective constant (mu) is P(x) = 0, where P(x) is:")
    print(f"({c4}) * mu^4 + ({c3}) * mu^3 + ({c2}) * mu^2 + ({c1}) * mu + ({c0}) = 0")
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # The connective constant is the unique positive real root
    connective_constant = None
    for root in roots:
        # Check if the root is real and positive
        if np.isreal(root) and root > 0:
            connective_constant = np.real(root)
            break
            
    print("\nNumerically, the roots of this polynomial are:")
    for r in roots:
        print(f"{r:.6f}")

    if connective_constant is not None:
        print(f"\nThe connective constant is the unique positive real root: mu = {connective_constant:.10f}")
    else:
        print("\nCould not find a positive real root for the connective constant.")

solve_connective_constant_polynomial()