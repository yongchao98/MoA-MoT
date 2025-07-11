import numpy as np
from scipy.linalg import svd

def construct_mandelbrot_matrix(n):
    """
    Constructs the Mandelbrot Matrix M_n of size (2^(n+1)-1) x (2^(n+1)-1).
    """
    N = 2**(n + 1) - 1
    M = np.zeros((N, N))

    # Superdiagonal of 1s
    if N > 1:
        np.fill_diagonal(M[:, 1:], 1)

    # Subdiagonal elements based on the Mandelbrot recurrence
    for i in range(1, N + 1):  # Using 1-based indexing for clarity
        if i % 2 == 0:  # i = 2k
            k = i // 2
            j = k
            M[i - 1, j - 1] = 1
        elif i > 1:  # i = 2k - 1
            k = (i + 1) // 2
            if k > 1:
                j = k
                M[i - 1, j - 1] = -1
    
    # Final element M_N,N = -1
    if N > 0:
        M[N - 1, N - 1] = -1
        
    return M

def ldlt_decomposition_diag(A):
    """
    Performs LDL' decomposition and returns the diagonal elements of D.
    This version is for symmetric positive-definite matrices.
    """
    N = A.shape[0]
    D_diag = np.zeros(N)
    L = np.eye(N)
    
    for j in range(N):
        v = np.zeros(j)
        for k in range(j):
            v[k] = L[j, k] * D_diag[k]
        
        D_diag[j] = A[j, j] - L[j, :j] @ v
        
        if abs(D_diag[j]) < 1e-12:
            return None # Matrix is singular or decomposition is unstable

        d_inv = 1.0 / D_diag[j]
        for i in range(j + 1, N):
            L[i, j] = (A[i, j] - L[i, :j] @ v) * d_inv
            
    return D_diag

def find_n0():
    """
    Finds the integer n_0 that minimizes Tr(D_n) * (Det(D_n))^(1/n).
    """
    print("Step 1: Finding the value of n0")
    min_f_val = float('inf')
    n0 = -1
    
    for n_test in range(1, 5):
        N = 2**(n_test + 1) - 1
        print(f"Testing n = {n_test} (matrix size: {N}x{N})")
        M_n = construct_mandelbrot_matrix(n_test)
        
        # We interpret "symmetric part" as M * M^T to ensure positive-definiteness
        S_n = M_n @ M_n.T
        
        D_n_diag = ldlt_decomposition_diag(S_n)
        
        if D_n_diag is None:
            print(f"  - LDL' decomposition failed for n={n_test}.")
            continue
            
        tr_D = np.sum(D_n_diag)
        det_D = np.prod(D_n_diag)
        
        # With S_n = M_n * M_n^T, det(S_n) = det(M_n)^2 >= 0, so this is safe.
        f_n = tr_D * (det_D**(1.0/n_test))
        
        print(f"  - Tr(D_{n_test}) = {tr_D:.4f}")
        print(f"  - Det(D_{n_test}) = {det_D:.4f}")
        print(f"  - F({n_test}) = {f_n:.4f}")

        if f_n < min_f_val:
            min_f_val = f_n
            n0 = n_test
            
    print(f"\nMinimum value of F(n) found at n0 = {n0}\n")
    return n0

def main():
    """
    Main function to solve the problem.
    """
    # Part 1: Find n0
    n0 = find_n0()

    if n0 == -1:
        print("Could not determine n0. Aborting.")
        return

    # Part 2: Calculations for n0
    print(f"Step 2: Performing calculations for n0 = {n0}")

    # Cofactor matrix of M_n0
    M = construct_mandelbrot_matrix(n0)
    print(f"M_{n0} is:\n{M}")
    
    det_M = np.linalg.det(M)
    if abs(det_M) < 1e-9:
        print("M_n0 is singular, cofactor matrix calculation using inverse is not possible this way.")
        return
        
    Cofactor_M = (np.linalg.inv(M) * det_M).T
    print(f"\nCofactor matrix C_{n0}:\n{Cofactor_M}")
    
    # Antisymmetric part of the cofactor matrix
    A_C = 0.5 * (Cofactor_M - Cofactor_M.T)
    print(f"\nAntisymmetric part A_C:\n{A_C}")

    # Permute A_C to get a tridiagonal matrix T
    # For n0=1, the matrix is 3x3. We can permute row/col 1 and 2 (0-indexed)
    P = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    T = P @ A_C @ P.T
    print(f"\nTridiagonal matrix T (by permuting A_C):\n{T}")

    # Square of T
    T_sq = T @ T
    print(f"\nSquare of the tridiagonal matrix, T^2:\n{T_sq}")
    
    # Singular values of T^2
    singular_values = svd(T_sq, compute_uv=False)
    print(f"\nSingular values of T^2 are: {singular_values}")
    
    # Largest Ky Fan norm (interpreted as the nuclear norm)
    nuclear_norm = np.sum(singular_values)
    
    print("\nStep 3: Calculating the final answer")
    print("The largest Ky Fan norm is the sum of all singular values (nuclear norm).")
    
    equation_str = " + ".join([f"{s:.4f}" for s in singular_values])
    print(f"Final calculation: {equation_str} = {nuclear_norm:.4f}")
    
    # Final answer in the required format
    print(f"\n<<< {nuclear_norm} >>>")

if __name__ == "__main__":
    main()
