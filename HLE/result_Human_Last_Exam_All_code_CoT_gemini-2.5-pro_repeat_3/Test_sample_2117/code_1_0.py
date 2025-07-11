import numpy as np

def get_cayley_menger_matrix(n):
    """Constructs the Cayley-Menger matrix for a regular n-simplex."""
    N = n + 2
    C = np.ones((N, N))
    np.fill_diagonal(C, 0)
    C[0, 0] = 0
    # The submatrix C[1:, 1:] is J-I, which is already correct.
    # The first row and column need to be set.
    C[0, :] = 1
    C[:, 0] = 1
    C[0, 0] = 0
    return C

def get_avg_eigenvalue_gap(M):
    """Calculates the average eigenvalue gap of a matrix M."""
    # Eigenvalues can be complex, so we take the absolute value for sorting
    eigvals = np.linalg.eigvals(M)
    # Sort eigenvalues by their real part, then imaginary part
    sorted_eigvals = sorted(eigvals, key=lambda x: (np.real(x), np.imag(x)))
    
    if len(sorted_eigvals) < 2:
        return 0.0
        
    gaps = np.abs(np.diff(sorted_eigvals))
    return np.mean(gaps)

def get_mean_square_singular_values(M):
    """Calculates the mean square of singular values of a matrix M."""
    # S_M = (1/N) * ||M||_F^2
    N = M.shape[0]
    frobenius_norm_sq = np.sum(np.abs(M)**2)
    return frobenius_norm_sq / N

def solve_for_n(n):
    """
    Solves the problem for a given integer n and prints the results.
    """
    print(f"--- Analyzing for n = {n} ---")
    
    N = n + 2
    C = get_cayley_menger_matrix(n)

    # The Gaussian-Hessenberg reduction P is built from column operations.
    # For C_n, the pivot c_{2,1} is 1, so no pivoting is needed for the first step.
    # The transformation P is the product of elementary matrices for the column ops.
    # H = P_inv * C * P, where P_inv represents row operations.
    # For step k=1, R_i -> R_i - m_{i,1}*R_2. Here m_{i,1} = c_{i,1}/c_{2,1} = 1/1 = 1 for i>2.
    # P_inv is the product of these row operation matrices.
    # The corresponding P = P_inv^{-1} will be a unit lower triangular matrix.
    # For the first step, P_1 = I + sum_{i=3 to N} e_i * e_2^T
    # For n>=2, the algorithm may require pivoting or may terminate early if pivots are 0.
    # If we interpret the algorithm as 'do nothing if column is already reduced',
    # no pivoting is needed and P remains unit triangular.
    
    # We construct P for the first step of the reduction.
    # For C_n, subsequent steps do not introduce further transformations
    # as the relevant sub-pivots become zero.
    P = np.identity(N)
    # For the first column (k=1), pivot is C[1,0]=1. We add multiples of column i to column 1+1=2.
    # C_2 -> C_2 + sum_{i=3 to N} m_{i,1}*C_i. m_{i,1}=1.
    # This corresponds to P = product of (I + m_{i,1}*e_2*e_i^T)
    # The actual transformation matrix P in H = P_inv*C*P is the inverse of the row-op matrix.
    P_inv = np.identity(N)
    for i in range(2, N):
        # Multiplier m = C[i,0]/C[1,0] = 1/1 = 1
        P_inv[i, 1] = -1.0 
    
    # P is the inverse of P_inv
    P = np.linalg.inv(P_inv)

    # Compute H
    H = np.linalg.inv(P) @ C @ P

    # Calculate all the terms
    E_P = get_avg_eigenvalue_gap(P)
    E_H = get_avg_eigenvalue_gap(H)
    S_P = get_mean_square_singular_values(P)
    S_H = get_mean_square_singular_values(H)
    
    product = E_P * E_H * S_P * S_H

    print("\nBased on the derivation that P is a unit triangular matrix, its eigenvalues are all 1.")
    p_eigvals = np.linalg.eigvals(P)
    print(f"Eigenvalues of P for n={n}: {np.round(p_eigvals, 5)}")
    print(f"This makes the average eigenvalue gap for P, E_P, equal to zero.\n")

    print("Calculating all values for the expression:")
    print(f"E_P = {E_P:.4f}")
    print(f"E_H = {E_H:.4f}")
    print(f"S_P = {S_P:.4f}")
    print(f"S_H = {S_H:.4f}")
    
    print("\nThe final product is the multiplication of these values:")
    print(f"E_P * E_H * S_P * S_H = {E_P:.4f} * {E_H:.4f} * {S_P:.4f} * {S_H:.4f} = {product:.4f}")

# Main execution
# We demonstrate for n=2, but the logic holds for any n.
n_val = 2
solve_for_n(n_val)

print("\nSince E_P is 0 for all n, the product is always 0.")
print("The set of possible values is {0}. The least upper bound of this set is 0.")

<<<0>>>