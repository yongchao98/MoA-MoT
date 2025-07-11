import numpy as np

def solve():
    """
    Calculates the infimum and supremum of the given expression.
    """
    # Define the polynomial f(x) = x^3 + x^2 + x - 1
    poly_coeffs = [1, 1, 1, -1]

    # 1. Find the roots of the polynomial
    roots = np.roots(poly_coeffs)

    # Identify the real root xi and one of the complex roots xi_1
    real_mask = np.isreal(roots)
    xi = roots[real_mask][0].real
    
    complex_mask = np.iscomplex(roots)
    xi_1 = roots[complex_mask][0]

    # 2. Define the eigenvectors u_1 and w_1
    # Eigenvector u_1 of M for eigenvalue xi_1
    # M = [[0, 0, 1], [1, 0, -1], [0, 1, -1]]
    # From M*u = xi_1*u, we can derive the eigenvector up to a scaling factor.
    # Let u_1 = [x, y, z]. We can set x=1.
    # z = xi_1 * x = xi_1
    # y = (x - z) / xi_1 = (1 - xi_1) / xi_1
    u_1 = np.array([1, (1 - xi_1) / xi_1, xi_1], dtype=np.complex128)

    # Eigenvector w_1 of M.T for eigenvalue xi_1
    # M.T is the companion matrix of f(x). Its eigenvector is known.
    w_1 = np.array([1, xi_1, xi_1**2], dtype=np.complex128)

    # 3. Calculate the terms for the asymptotic formula
    # |w_1^T u_1|^2
    w1_dot_u1 = np.dot(w_1, u_1)
    abs_sq_w1_dot_u1 = np.abs(w1_dot_u1)**2

    # ||u_1||^2
    u1_norm_sq = np.linalg.norm(u_1)**2

    # |u_1 . u_1|
    u1_dot_u1 = np.dot(u_1, u_1)
    abs_u1_dot_u1 = np.abs(u1_dot_u1)

    # 4. Calculate the infimum and supremum
    inf_val = (2 * u1_norm_sq - 2 * abs_u1_dot_u1) / abs_sq_w1_dot_u1
    sup_val = (2 * u1_norm_sq + 2 * abs_u1_dot_u1) / abs_sq_w1_dot_u1
    
    print(f"The infimum is: {inf_val}")
    print(f"The supremum is: {sup_val}")

solve()