import numpy as np

def solve():
    """
    Solves the problem by finding the eigenvalues of the associated quadratic form matrix Q.
    """
    # Step 1: Find the real root xi of the polynomial f(x) = x^3 + x^2 + x - 1.
    # The coefficients are given in descending order of power.
    poly_f_coeffs = [1, 1, 1, -1]
    roots_f = np.roots(poly_f_coeffs)
    
    # The real root is needed for calculations.
    xi = None
    for r in roots_f:
        if np.isreal(r):
            xi = r.real
            break

    # Step 2: Construct the matrix Q = V^H V.
    # The elements of Q can be calculated from sums of powers of the roots,
    # which can be simplified to expressions in terms of xi.
    # Q_ij = sum_{k} conj(xi_k^{i-1}) * xi_k^{j-1}
    Q = np.zeros((3, 3))
    
    # Q_11 = sum |1|^2 = 3
    Q[0, 0] = 3.0
    
    # Q_12 = sum xi_k = -1 (from Vieta's formulas)
    # Q_21 = sum conj(xi_k) = -1
    Q[0, 1] = Q[1, 0] = -1.0
    
    # Q_13 = sum xi_k^2 = (sum xi_k)^2 - 2*sum(xi_i*xi_j) = (-1)^2 - 2*(1) = -1
    # Q_31 = sum conj(xi_k^2) = -1
    Q[0, 2] = Q[2, 0] = -1.0
    
    # Q_22 = sum |xi_k|^2 = xi^2 + |xi_1|^2 + |xi_2|^2 = xi^2 + 2/xi
    Q[1, 1] = xi**2 + 2 / xi
    
    # Q_23 = sum conj(xi_k)*xi_k^2 = sum |xi_k|^2 * xi_k
    #      = xi^3 + (1/xi)*xi_1 + (1/xi)*xi_2 = xi^3 + (1/xi)*(xi_1+xi_2)
    #      = xi^3 + (1/xi)*(-1-xi) = (1-xi-xi^2) - 1/xi - 1 = -xi - xi^2 - 1/xi
    # Since Q is real, Q_32 = Q_23
    Q[1, 2] = Q[2, 1] = -xi - xi**2 - 1 / xi
    
    # Q_33 = sum |xi_k^2|^2 = xi^4 + |xi_1|^4 + |xi_2|^4 = xi^4 + 2/xi^2
    Q[2, 2] = xi**4 + 2 / xi**2

    # Step 3: Compute the eigenvalues of Q.
    # Since Q is real and symmetric, its eigenvalues are real.
    eigenvalues = np.linalg.eigvalsh(Q)
    lambda_min = np.min(eigenvalues)
    lambda_max = np.max(eigenvalues)

    # Step 4: The asymptotic bounds for E_n are 2/lambda_max and 2/lambda_min.
    # These bounds give the infimum and supremum over all n.
    inf_val = 2 / lambda_max
    sup_val = 2 / lambda_min

    # Print the results as requested.
    print(f"inf_n|P_n(xi)|(a_n^2+b_n^2+c_n^2) = {inf_val}")
    print(f"sup_n|P_n(xi)|(a_n^2+b_n^2+c_n^2) = {sup_val}")

solve()