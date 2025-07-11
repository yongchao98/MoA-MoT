import numpy as np

def solve_polynomial_optimization():
    """
    Solves the problem of finding the infimum and supremum of the given expression.
    
    The problem considers a polynomial f(x) = x^3 + x^2 + x - 1, with roots xi (real)
    and xi_1, xi_2 (complex). For any positive integer n, we have a polynomial
    P_n(x) = a_n + b_n*x + c_n*x^2 such that P_n(root) = root^n for each root.
    The goal is to find the infimum and supremum of the expression:
    |P_n(xi)| * (a_n^2 + b_n^2 + c_n^2).
    
    Through analytical derivation, the supremum and infimum are found to be functions
    of the real root xi.
    """
    
    # Coefficients of the polynomial x^3 + x^2 + x - 1
    coeffs = [1, 1, 1, -1]
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # Identify the real root xi
    xi = None
    for root in roots:
        if np.isreal(root):
            xi = np.real(root)
            break
            
    if xi is None:
        print("Could not find the real root.")
        return
        
    xi_sq = xi**2
    
    # Formulas for infimum and supremum derived analytically
    # Supremum = (7 + 4*xi + 3*xi^2) / 11
    # Infimum = (1 + 8*xi + 5*xi^2) / 11
    
    sup_val = (7 + 4 * xi + 3 * xi_sq) / 11
    inf_val = (1 + 8 * xi + 5 * xi_sq) / 11
    
    print("The real root is xi ≈ {:.6f}".format(xi))
    print("xi^2 ≈ {:.6f}".format(xi_sq))
    print("\n--- Finding the infimum ---")
    print("inf_n |P_n(xi)|(a_n^2+b_n^2+c_n^2) = (1 + 8*xi + 5*xi^2) / 11")
    print("= (1 + 8*{:.6f} + 5*{:.6f}) / 11".format(xi, xi_sq))
    print("= {:.6f}".format(inf_val))
    
    print("\n--- Finding the supremum ---")
    print("sup_n |P_n(xi)|(a_n^2+b_n^2+c_n^2) = (7 + 4*xi + 3*xi^2) / 11")
    print("= (7 + 4*{:.6f} + 3*{:.6f}) / 11".format(xi, xi_sq))
    print("= {:.6f}".format(sup_val))

solve_polynomial_optimization()