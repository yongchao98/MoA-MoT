import numpy as np

def solve_polynomial_problem():
    """
    Solves the problem by calculating initial values and asymptotic limits
    of the sequence E_n.
    """
    # 1. Find the roots of the polynomial f(x) = x^3 + x^2 + x - 1
    poly_coeffs = [1, 1, 1, -1]
    roots = np.roots(poly_coeffs)
    
    xi = roots[np.isreal(roots)][0].real
    complex_roots = roots[np.iscomplex(roots)]
    xi1 = complex_roots[0]
    xi2 = np.conj(xi1)

    # 2. Asymptotic analysis for large n
    # The limit points of E_n are determined by the eigenvectors and eigenvalues.
    # The vector v_n can be expressed as v_n = K1*xi^n + K2*xi1^n + K3*xi2^n,
    # where Ki are columns of the inverse Vandermonde matrix.
    V = np.array([[1, xi, xi**2], [1, xi1, xi1**2], [1, xi2, xi2**2]])
    try:
        V_inv = np.linalg.inv(V)
    except np.linalg.LinAlgError:
        print("Vandermonde matrix is singular.")
        return

    # K2 is the coefficient vector for the dominant term xi1^n
    K2 = V_inv[:, 1]
    
    # For large n, E_n approaches an oscillating value.
    # The mean of this oscillation is 2 * ||K2||^2
    # The amplitude is 2 * |K2.T @ K2|
    norm_K2_sq = np.linalg.norm(K2)**2
    K2_T_K2_abs = np.abs(K2.T @ K2)
    
    lim_sup = 2 * norm_K2_sq + 2 * K2_T_K2_abs
    lim_inf = 2 * norm_K2_sq - 2 * K2_T_K2_abs

    # 3. Calculate initial values of E_n for small n
    # We use the recurrence v_{n+1} = M @ v_n
    # M is the companion matrix for f(x) = x^3 + x^2 + x - 1
    M = np.array([[0, 0, 1], [1, 0, -1], [0, 1, -1]], dtype=float)
    
    # Initial vector for n=1 (P_1(x) = x)
    v = np.array([0, 1, 0], dtype=float)
    
    num_initial_terms = 20
    e_values = []
    
    for n in range(1, num_initial_terms + 1):
        norm_v_sq = np.linalg.norm(v)**2
        e_n = (xi**n) * norm_v_sq
        e_values.append(e_n)
        # Update v for the next iteration
        v = M @ v

    # 4. Determine the overall infimum and supremum
    overall_sup = max(max(e_values), lim_sup)
    overall_inf = min(min(e_values), lim_inf)
    
    print(f"The infimum of |P_n(xi)|(a_n^2+b_n^2+c_n^2) is: {overall_inf}")
    print(f"The supremum of |P_n(xi)|(a_n^2+b_n^2+c_n^2) is: {overall_sup}")

solve_polynomial_problem()