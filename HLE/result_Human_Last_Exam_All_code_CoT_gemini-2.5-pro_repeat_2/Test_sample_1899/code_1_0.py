import numpy as np

def solve():
    """
    Solves the problem by numerically calculating the exact theoretical
    expressions for the infimum and supremum.
    """
    # Define the polynomial coefficients for f(x) = x^3 + x^2 + x - 1
    coeffs = [1, 1, 1, -1]
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # Assign roots to xi, xi_1, xi_2
    # The real root is xi. The problem states xi is approx 0.54.
    real_root_mask = np.isreal(roots)
    xi = roots[real_root_mask][0].real
    
    complex_roots = roots[~real_root_mask]
    # Order of complex roots doesn't matter as they are conjugates
    xi_1 = complex_roots[0]
    xi_2 = complex_roots[1]

    # f'(x) = 3x^2 + 2x + 1
    def f_prime(x):
        return 3*x**2 + 2*x + 1

    # Calculate the necessary components for the formula
    
    # |f'(xi_1)|^2
    f_prime_xi_1 = f_prime(xi_1)
    abs_f_prime_xi_1_sq = np.abs(f_prime_xi_1)**2
    
    # |1 + xi_1|^2
    one_plus_xi_1 = 1 + xi_1
    abs_one_plus_xi_1_sq = np.abs(one_plus_xi_1)**2
    
    # The complex number Z_0 = 1 + xi^2 * xi_2^2 + (1 + xi_1)^2
    Z_0 = 1 + (xi**2) * (xi_2**2) + (1 + xi_1)**2
    abs_Z_0 = np.abs(Z_0)
    
    # The constant part of the expression
    C1 = xi + abs_one_plus_xi_1_sq + 1
    
    # The scaling factor
    factor = 2 / abs_f_prime_xi_1_sq

    # Calculate infimum and supremum
    inf_value = factor * (C1 - abs_Z_0)
    sup_value = factor * (C1 + abs_Z_0)

    print(f"The infimum is: {inf_value}")
    print(f"The supremum is: {sup_value}")

solve()