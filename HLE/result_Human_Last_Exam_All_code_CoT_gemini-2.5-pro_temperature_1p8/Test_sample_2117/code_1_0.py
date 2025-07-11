import numpy as np
from scipy.linalg import hessenberg

def calculate_metrics(M):
    """Calculates the average eigenvalue gap (E_M) and mean square of singular values (S_M) for a matrix M."""
    k = M.shape[0]

    # Eigenvalues
    eigvals = np.linalg.eigvals(M)
    # For E_M, we use the spread of the real parts of the eigenvalues
    real_eigvals = np.real(eigvals)
    E_M = (np.max(real_eigvals) - np.min(real_eigvals)) / (k - 1) if k > 1 else 0

    # Singular values
    sigmas_sq = np.linalg.svd(M, compute_uv=False)**2
    S_M = np.mean(sigmas_sq)

    return E_M, S_M

def solve_for_n(n):
    """
    Constructs the Cayley-Menger matrix for a given n, decomposes it,
    and calculates the product E_P*E_H*S_P*S_H.
    """
    k = n + 2

    # 1. Construct the Cayley-Menger matrix C_n
    if n == 0: # A 0-simplex is a single point. C is 2x2.
        C = np.array([[0., 1.], [1., 0.]])
    else:
        C = np.zeros((k, k))
        # Top-left block for vertices
        J_minus_I = np.ones((n + 1, n + 1)) - np.eye(n + 1)
        C[1:, 1:] = J_minus_I
        # Border of 1s
        C[0, 1:] = 1.0
        C[1:, 0] = 1.0

    # 2. Perform Hessenberg Decomposition: C = P H P^T
    # scipy.linalg.hessenberg returns H, P such that H = P.T @ C @ P
    # So C = P @ H @ P.T. P is orthogonal.
    H, P = hessenberg(C, calc_q=True)
    
    # 3. Calculate metrics for P and H
    E_P, S_P = calculate_metrics(P)
    E_H, S_H = calculate_metrics(H)
    
    # 4. Compute the product
    product = E_P * S_P * E_H * S_H

    # Also compute theoretical values for comparison
    E_H_theory = (n + 2) / (n + 1) if n > 0 else (0-(-2))/(2-1) # special case n=0 not requested
    S_H_theory = float(n + 1)
    
    # Print the equation with calculated numbers
    print(f"For n = {n}:")
    print(f"  E_P = {E_P:.4f}, S_P = {S_P:.4f}")
    print(f"  E_H = {E_H:.4f} (Theory: {E_H_theory:.4f})")
    print(f"  S_H = {S_H:.4f} (Theory: {S_H_theory:.4f})")
    print(f"  Product E_P*E_H*S_P*S_H = {product:.6f}")
    print(f"  Theoretical Bound (2 + 2/(n+1)) = {2 + 2/(n+1):.6f}")
    print("-" * 20)


if __name__ == "__main__":
    # Calculate for a range of n values
    for n_val in range(1, 11):
        solve_for_n(n_val)
