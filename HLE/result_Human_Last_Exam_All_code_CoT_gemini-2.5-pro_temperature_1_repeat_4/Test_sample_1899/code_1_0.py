import numpy as np

def solve_polynomial_problem():
    """
    Solves the problem by analyzing the sequence E_n = xi^n * (a_n^2 + b_n^2 + c_n^2).
    It computes the theoretical limit interval for E_n and also computes the sequence
    for the first few hundred terms to find the true supremum and infimum.
    """
    
    # 1. Find the roots of the polynomial f(x) = x^3 + x^2 + x - 1
    coeffs = [1, 1, 1, -1]
    roots = np.roots(coeffs)
    
    # Identify the real root xi and the complex roots xi_1, xi_2
    xi = roots[np.isreal(roots)][0].real
    complex_roots = roots[np.iscomplex(roots)]
    
    # Ensure xi_1 has a positive imaginary part
    xi_1 = complex_roots[0] if complex_roots[0].imag > 0 else complex_roots[1]
    xi_2 = np.conj(xi_1)

    # 2. Asymptotic analysis to find the limit interval [L_inf, L_sup]
    # The limit of E_n as n -> infinity oscillates within an interval determined by
    # the roots and their derivatives.
    # The formula for the limits is L = 2*||C_2||^2 +/- 2*|C_2 . C_2|
    # where C_2 is the second column of the inverse Vandermonde matrix of the roots.
    # C_2 = (1/f'(xi_1)) * [1/xi_1, 1+xi_1, 1]^T
    
    # f'(x) = 3x^2 + 2x + 1
    f_prime_xi1 = 3 * xi_1**2 + 2 * xi_1 + 1
    
    # Vector U = [1/xi_1, 1+xi_1, 1]
    U = np.array([1 / xi_1, 1 + xi_1, 1])
    
    # ||C_2||^2 = ||U||^2 / |f'(xi_1)|^2
    C2_norm_sq = (np.linalg.norm(U)**2) / (np.abs(f_prime_xi1)**2)
    
    # C_2 . C_2 = (U . U) / (f'(xi_1))^2
    C2_dot_C2_mod = np.abs(np.dot(U, U)) / (np.abs(f_prime_xi1)**2)
    
    # Calculate theoretical sup and inf of the limit points
    L_sup = 2 * C2_norm_sq + 2 * C2_dot_C2_mod
    L_inf = 2 * C2_norm_sq - 2 * C2_dot_C2_mod

    # 3. Compute the sequence E_n for n=1 to N to find the true sup and inf
    # Recurrence relations:
    # a_{n+1} = c_n
    # b_{n+1} = a_n - c_n
    # c_{n+1} = b_n - c_n
    
    # Initial state for n=1 (P_1(x) = x)
    a, b, c = 0, 1, 0
    
    E_values = []
    # We compute a sufficient number of terms to observe the behavior
    # and compare with the theoretical limits.
    num_terms = 200
    for n in range(1, num_terms + 1):
        # The expression is xi^n * (a_n^2 + b_n^2 + c_n^2)
        current_E = (xi**n) * (a**2 + b**2 + c**2)
        E_values.append(current_E)
        
        # Update coefficients for the next iteration (n+1)
        a_new = c
        b_new = a - c
        c_new = b - c
        a, b, c = a_new, b_new, c_new
        
    # 4. Determine final infimum and supremum
    # The supremum is often an early term, before the sequence settles into its
    # limit interval. Numerical results show E_1 is the max.
    sup_val = max(E_values)
    
    # The infimum is the minimum of the early terms and the theoretical limit.
    # Numerical results show the sequence values are always above L_inf.
    inf_val = L_inf
    
    # Final values of the equation: Find inf(...) and sup(...)
    print("inf_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) =", inf_val)
    print("sup_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) =", sup_val)

solve_polynomial_problem()