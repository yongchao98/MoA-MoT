import numpy as np
from scipy.linalg import hessenberg

def solve():
    """
    This function calculates the product for n=1, which corresponds to the
    least upper bound of the expression over all positive integers n.
    """
    n = 1
    m = n + 2  # Size of the matrix is (n+2) x (n+2)

    # 1. Construct the Cayley-Menger matrix C_n for a regular n-simplex
    # For n=1, this is a 3x3 matrix.
    C_n = np.ones((m, m))
    C_n[1:, 1:] -= np.identity(n + 1)
    C_n[0, 0] = 0

    # 2. Perform Hessenberg decomposition: C_n = P * H * P^T
    # We assume an orthogonal decomposition, where P is orthogonal.
    # scipy.linalg.hessenberg returns H and P (as Q) such that H = P^T * C_n * P
    # So C_n = P * H * P^T, which matches our formulation.
    H, P = hessenberg(C_n, calc_q=True)

    # 3. Calculate E_H and S_H
    # Eigenvalues of H are the same as C_n. C_n is symmetric, so we use eigvalsh.
    eigvals_H = np.sort(np.linalg.eigvalsh(C_n))
    
    # E_H: Average eigenvalue gap of H
    E_H = (eigvals_H[-1] - eigvals_H[0]) / (m - 1)
    
    # S_H: Mean square of singular values of H. For a symmetric matrix,
    # singular values are the absolute values of eigenvalues.
    s_vals_sq_H = eigvals_H**2
    S_H = np.mean(s_vals_sq_H)

    # 4. Calculate E_P and S_P
    # P is an orthogonal matrix.
    # S_P: Mean square of singular values of P. For an orthogonal matrix, all are 1.
    S_P = 1.0
    
    # E_P: Average eigenvalue gap of P.
    eigvals_P = np.linalg.eigvals(P)
    real_parts_P = np.sort(np.real(eigvals_P))
    E_P = (real_parts_P[-1] - real_parts_P[0]) / (m - 1)
    
    # 5. Compute the final product
    product = E_P * E_H * S_P * S_H

    # Print the final equation with the computed values
    print(f"For n = {n}, the components are:")
    print(f"E_P = {E_P:.4f}")
    print(f"S_P = {S_P:.4f}")
    print(f"E_H = {E_H:.4f}")
    print(f"S_H = {S_H:.4f}")
    print(f"The product is {E_P:.4f} * {S_P:.4f} * {E_H:.4f} * {S_H:.4f} = {product:.4f}")

solve()