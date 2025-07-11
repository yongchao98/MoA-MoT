import numpy as np

def solve():
    """
    Solves the problem by finding the infimum and supremum of the given expression.
    """
    # 1. Find the roots of the polynomial f(x) = x^3 + x^2 + x - 1
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    
    # Extract the real root xi
    xi = None
    for r in roots:
        if np.isreal(r):
            xi = r.real
            break
            
    if xi is None:
        print("Could not find the real root.")
        return

    # 2. Calculate the supremum
    # Based on analysis of the first few terms, the supremum is E_1 = xi.
    # a_1=0, b_1=1, c_1=0.
    # E_1 = xi^1 * (0^2 + 1^2 + 0^2) = xi
    sup_value = xi
    
    # 3. Calculate the infimum using the asymptotic formula
    # inf_value = (3*xi^2 + 2*xi + 1) * (1 - xi)^2 / 22
    term1 = 3 * xi**2 + 2 * xi + 1
    term2 = (1 - xi)**2
    inf_value = (term1 * term2) / 22
    
    # Print the results
    print(f"The root xi is approximately: {xi}")
    print(f"The infimum value is: {inf_value}")
    print(f"The supremum value is: {sup_value}")

solve()
