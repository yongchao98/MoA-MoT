import numpy as np

def solve_connective_constant_poly():
    """
    This function provides the minimal polynomial for the connective constant
    of the specified graph G and calculates the constant's value.
    
    The minimal polynomial is known from literature on lattice combinatorics:
    P(x) = x^3 - 2x^2 - 2x - 2 = 0
    The connective constant mu is the largest real root of this polynomial.
    """
    
    # Coefficients of the polynomial P(x) = c3*x^3 + c2*x^2 + c1*x + c0
    c3 = 1
    c2 = -2
    c1 = -2
    c0 = -2
    
    coeffs = [c3, c2, c1, c0]
    
    # We print the minimal polynomial equation as requested, showing each coefficient.
    print("The minimal polynomial for the connective constant of graph G is P(x) = 0, where:")
    print(f"P(x) = ({c3})*x^3 + ({c2})*x^2 + ({c1})*x + ({c0})")
    
    # The Rational Root Theorem can be used to show this polynomial is irreducible
    # over Q, confirming it is the minimal polynomial.
    
    # For completeness, we calculate the numerical value of the connective constant.
    try:
        roots = np.roots(coeffs)
        
        # The connective constant mu is the single positive real root.
        mu = 0.0
        for r in roots:
            if np.isreal(r):
                real_root = np.real(r)
                if real_root > mu:
                    mu = real_root
                    
        print(f"\nThe connective constant mu is the largest real root of P(x), which is approximately: {mu:.9f}")
        
    except ImportError:
        print("\nNote: 'numpy' library is not installed. Cannot compute the numerical value of the root.")
        print("To install it, run: pip install numpy")


if __name__ == "__main__":
    solve_connective_constant_poly()