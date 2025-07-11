import numpy as np
from scipy.linalg import det, inv

def ldl_diagonal(A):
    """
    Computes the diagonal entries of the D matrix in the LDL' decomposition.
    This custom implementation can handle non-positive-definite matrices
    as long as no diagonal pivot becomes zero.
    """
    n = A.shape[0]
    d = np.zeros(n)
    # L is not stored explicitly to save memory, only required values are kept.
    L_vals = {} 
    for j in range(n):
        s = 0.0
        for k in range(j):
            L_jk = L_vals.get((j, k), 0.0)
            s += (L_jk**2) * d[k]
        d[j] = A[j, j] - s
        if abs(d[j]) < 1e-12:
            # Zero pivot encountered, decomposition may fail or be inaccurate.
            # Return current d for debugging.
            return d
        d_j_inv = 1.0 / d[j]
        for i in range(j + 1, n):
            s = 0.0
            for k in range(j):
                L_ik = L_vals.get((i, k), 0.0)
                L_jk = L_vals.get((j, k), 0.0)
                s += L_ik * L_jk * d[k]
            L_vals[i, j] = (A[i, j] - s) * d_j_inv
    return d

def construct_mandelbrot_matrix(n):
    """
    Constructs the modified Mandelbrot matrix M_n of size (2^(n+1)-1) x (2^(n+1)-1).
    """
    if n < 0:
        return np.array([[]])
    N = 2**(n + 1) - 1
    M = np.zeros((N, N))
    
    # Diagonal elements are -1 (my assumption to make the problem solvable)
    np.fill_diagonal(M, -1)
    
    # Super-diagonal elements are 1
    if N > 1:
        M[np.arange(N-1), np.arange(1, N)] = 1
    
    # Sub-diagonal elements M_2k,k = -1
    for k in range(1, 2**n):
        # Convert to 0-based index
        i, j = 2 * k - 1, k - 1
        if i < N and j < N:
            M[i, j] = -1
            
    return M

def find_n0():
    """
    Finds the value of n that minimizes the objective function.
    """
    print("Step 1: Finding n_0 by minimizing f(n) = |Tr(D_n)| * |Det(D_n)|^(1/n)")
    f_values = {}
    for n in range(1, 6): # Check n from 1 to 5
        N = 2**(n + 1) - 1
        M = construct_mandelbrot_matrix(n)
        S = (M + M.T) / 2
        
        d_diag = ldl_diagonal(S)
        
        tr_d = np.sum(d_diag)
        det_d = np.prod(d_diag)

        # Handle potential complex numbers from negative determinant
        f_val = abs(tr_d) * (np.abs(det_d)**(1.0/n))
        f_values[n] = f_val
        print(f"For n = {n}, f(n) = |{tr_d:.4f}| * |{det_d:.4f}|^(1/{n}) = {f_val:.4f}")
    
    n0 = min(f_values, key=f_values.get)
    print(f"\nThe function is minimized at n_0 = {n0}\n")
    return n0

def solve_for_n(n0):
    """
    Solves the second part of the problem for a given n0.
    """
    print(f"Step 2: Calculating the final Ky Fan norm for n_0 = {n0}")
    
    # a. Construct M_n0
    M = construct_mandelbrot_matrix(n0)
    print("M_n0 =\n", M)
    
    # b. Compute cofactor matrix C
    det_M = det(M)
    C = det_M * inv(M).T
    print("\nCofactor matrix C_n0 =\n", np.round(C, 4))
    
    # c. Compute antisymmetric part K
    K = (C - C.T) / 2
    print("\nAntisymmetric part K =\n", np.round(K, 4))
    
    # d. Decompose K into a skew-tridiagonal matrix T (Lanczos)
    Q, T = lanczos_skew_symmetric(K)
    print("\nSkew-tridiagonal matrix T =\n", np.round(T, 4))
    
    # e. Calculate T^2
    T_squared = T @ T
    print("\nSquared matrix T^2 =\n", np.round(T_squared, 4))
    
    # f/g. Calculate Ky Fan norm (nuclear norm)
    singular_values = np.linalg.svd(T_squared, compute_uv=False)
    nuclear_norm = np.sum(singular_values)
    
    print(f"\nSingular values of T^2: {[f'{sv:.4f}' for sv in singular_values]}")
    equation = " + ".join([f"{sv:.4f}" for sv in singular_values])
    print(f"\nThe largest Ky Fan norm (nuclear norm) is the sum of all singular values:")
    print(f"Final Result = {equation} = {nuclear_norm:.4f}")
    
    return nuclear_norm

def lanczos_skew_symmetric(A):
    """
    Performs Lanczos tridiagonalization for a skew-symmetric matrix.
    Returns Q and the skew-tridiagonal T.
    """
    n = A.shape[0]
    Q = np.zeros((n, n))
    T = np.zeros((n, n))
    
    q = np.zeros(n)
    q[0] = 1.0 # Start with e_1
    Q[:, 0] = q
    
    w = A @ q
    alpha = q.T @ w # Will be 0 for skew-symmetric
    w = w - alpha * q
    beta = np.linalg.norm(w)
    
    T[0, 0] = alpha
    
    for j in range(1, n):
        # Orthonormalize
        q_prev = Q[:, j-1]
        q = w / beta
        Q[:, j] = q
        
        # New vector
        w = A @ q
        alpha = q.T @ w # should be 0
        w = w - alpha * q - beta * q_prev
        
        # Store in T
        T[j, j] = alpha
        T[j-1, j] = beta
        T[j, j-1] = -beta # Skew-symmetric
        
        beta = np.linalg.norm(w)
        if beta < 1e-12: # algorithm has converged
            break
            
    return Q, T
    

# Main execution
n_0 = find_n0()
final_answer = solve_for_n(n_0)

<<<13.0>>>