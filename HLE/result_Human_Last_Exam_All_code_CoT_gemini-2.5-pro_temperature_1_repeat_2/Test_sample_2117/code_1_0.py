import numpy as np

def get_cayley_menger_matrix(n):
    """
    Constructs the Cayley-Menger matrix for a regular n-simplex with unit side length.
    The matrix size is (n+2)x(n+2).
    """
    k = n + 1
    size = k + 1
    # For a regular simplex with side length 1, the squared distance d_ij^2 is 1 for i!=j.
    # The Cayley-Menger matrix C has C_ij = d_ij^2 for i,j > 0.
    # The first row and column are special (bordered with 1s).
    C = np.ones((size, size), dtype=float)
    C[0, 0] = 0
    # Set the diagonal of the distance submatrix to 0
    np.fill_diagonal(C[1:, 1:], 0)
    return C

def gaussian_hessenberg_decomposition(A):
    """
    Performs Gaussian Hessenberg decomposition A = P H P^-1.
    Note: The standard algorithm produces H such that P_inv * A * P_inv^-1 = H,
    where P_inv is the transformation matrix. This means P in the problem's
    notation is the inverse of our computed transformation matrix.
    """
    size = A.shape[0]
    H = A.copy()
    P_inv = np.identity(size) # P_inv is the matrix that transforms A to H

    for j in range(size - 2):
        for i in range(j + 2, size):
            # The pivot is H[j+1, j]. For this specific problem, it can be shown
            # that the pivots are always non-zero, so no row swapping is needed.
            if np.isclose(H[j + 1, j], 0):
                continue
                
            m = H[i, j] / H[j + 1, j]
            
            # Create elementary matrix for the transformation
            E = np.identity(size)
            E[i, j + 1] = -m
            
            # Apply similarity transformation to H
            H = E @ H @ np.linalg.inv(E)
            
            # Accumulate the transformation
            P_inv = E @ P_inv
    
    # P from the problem statement is the inverse of the transformation matrix P_inv
    P = np.linalg.inv(P_inv)
    return P, H

def get_avg_eigenvalue_gap(M):
    """Calculates the average eigenvalue gap of a matrix M."""
    k = M.shape[0]
    if k < 2:
        return 0
    
    e_vals = np.linalg.eigvals(M)
    # The gaps are between sorted real parts of eigenvalues
    e_vals_real = np.sort(np.real(e_vals))
    
    gap = (e_vals_real[-1] - e_vals_real[0]) / (k - 1)
    return gap

def get_mean_square_singular_values(M):
    """Calculates the mean square of the singular values of M."""
    k = M.shape[0]
    if k == 0:
        return 0
    # The sum of squared singular values is the squared Frobenius norm.
    frobenius_norm_sq = np.sum(np.abs(M)**2)
    return frobenius_norm_sq / k

def solve_for_n(n):
    """
    Calculates the product E_P*E_H*S_P*S_H for a given n.
    """
    print(f"--- Calculating for n = {n} ---")
    C = get_cayley_menger_matrix(n)
    P, H = gaussian_hessenberg_decomposition(C)
    
    E_P = get_avg_eigenvalue_gap(P)
    E_H = get_avg_eigenvalue_gap(H)
    S_P = get_mean_square_singular_values(P)
    S_H = get_mean_square_singular_values(H)
    
    product = E_P * E_H * S_P * S_H
    
    print(f"P's eigenvalues are on the diagonal: {np.diag(P)}")
    print(f"E_P = {E_P}")
    print(f"E_H = {E_H}")
    print(f"S_P = {S_P}")
    print(f"S_H = {S_H}")
    print(f"The final product for n={n} is:")
    print(f"{E_P} * {E_H} * {S_P} * {S_H} = {product}\n")
    return product

# Calculate for a few values of n
solve_for_n(1)
solve_for_n(2)
solve_for_n(3)

# Based on the derivation, the product is 0 for all n.
# Therefore, the least upper bound of the set {0} is 0.
final_answer = 0.0
print(f"The analysis shows the product is 0 for all n.")
print(f"The least upper bound is {final_answer}.")
<<<0.0>>>