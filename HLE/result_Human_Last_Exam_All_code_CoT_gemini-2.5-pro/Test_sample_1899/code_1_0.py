import numpy as np

def solve():
    """
    This function calculates the infimum and supremum of the given expression.
    """
    # The polynomial is f(x) = x^3 + x^2 + x - 1
    coeffs = [1, 1, 1, -1]
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # Identify the real root xi. It's the root with a zero imaginary part.
    xi = 0.0
    for r in roots:
        if np.isreal(r):
            xi = np.real(r)
            break
            
    # The infimum of the expression is 0.
    inf_val = 0
    
    # The supremum is derived from the asymptotic behavior of the coefficients.
    # The analytical formula for the supremum is 4*(1+xi) / (7 + 11*xi^2).
    sup_val = 4 * (1 + xi) / (7 + 11 * xi**2)
    
    # Output the results
    print("The infimum is:")
    print(f"inf |P_n(xi)| (a_n^2+b_n^2+c_n^2) = {inf_val}")
    print("\nThe supremum is:")
    print(f"sup |P_n(xi)| (a_n^2+b_n^2+c_n^2) = {sup_val}")

solve()