import numpy as np

def get_h_matrix(k):
    """
    Constructs the matrix H_k of size (2^k-1) x (2^k-1) recursively.
    These are a known family of matrices whose eigenvalues lie on the
    boundary of the Mandelbrot set.
    """
    if k == 1:
        return np.array([[-1.0]])
    
    h_prev = get_h_matrix(k - 1)
    m = h_prev.shape[0]  # Size of the previous matrix, 2^(k-1)-1
    size = 2 * m + 1    # Size of the current matrix, 2^k-1
    
    h_curr = np.zeros((size, size))
    
    # Place H_{k-1} in the top-left corner
    h_curr[0:m, 0:m] = h_prev
    # Place identity matrix
    h_curr[0:m, m:(2 * m)] = np.eye(m)
    # Place the single '1'
    h_curr[m, 2 * m] = 1
    # Place the negative identity matrix
    h_curr[(m + 1):size, 0:m] = -np.eye(m)
    
    return h_curr

def get_mandelbrot_matrix(n):
    """
    The problem defines M_n as a matrix of size (2^(n+1)-1).
    This corresponds to H_{n+1} in our recursive construction.
    """
    return get_h_matrix(n + 1)

def ldl_decomposition(A):
    """
    Computes the LDL' decomposition of a symmetric matrix A without pivoting.
    Returns the diagonal of the matrix D.
    """
    n = A.shape[0]
    d_diag = np.zeros(n)
    L = np.eye(n)
    
    for i in range(n):
        sum_ld = sum(L[i, k]**2 * d_diag[k] for k in range(i))
        d_diag[i] = A[i, i] - sum_ld
        
        if abs(d_diag[i]) < 1e-12:
            # Matrix is singular or near-singular, subsequent calculations invalid/unstable.
            # Set to exactly zero for determinant calculation.
            d_diag[i] = 0.0
            continue # Continue to avoid division by zero
            
        for j in range(i + 1, n):
            sum_ldl = sum(L[j, k] * L[i, k] * d_diag[k] for k in range(i))
            L[j, i] = (A[j, i] - sum_ldl) / d_diag[i]
            
    return d_diag

# --- Step 1 & 2: Find n_0 ---
f_values = {}
n_range = range(1, 5) # Check for n=1,2,3,4

for n in n_range:
    M_n = get_mandelbrot_matrix(n)
    S_n = 0.5 * (M_n + M_n.T)
    
    D_n_diag = ldl_decomposition(S_n)
    
    det_D_n = np.prod(D_n_diag)
    tr_D_n = np.sum(D_n_diag)
    
    # Objective function is Tr(D_n) * (Det(D_n))^(1/n)
    # The minimum is expected to be over real values. If Det is negative
    # and n is even, the result is complex, which we treat as infinity.
    if det_D_n < 0 and n % 2 == 0:
        f_values[n] = float('inf')
    else:
        # This computes real roots for negative numbers if n is odd.
        val = np.sign(det_D_n) * (np.abs(det_D_n)**(1.0/n))
        f = tr_D_n * val
        f_values[n] = f

# Find n_0 corresponding to the minimum value of f(n)
n0 = min(f_values, key=f_values.get)

# --- Step 3: Process M_n0 ---
M = get_mandelbrot_matrix(n0)

# Compute cofactor matrix C = det(M) * (M^-1)^T
# Our matrix construction gives det(M) = -1 for all n>=1
det_M = -1.0 
inv_M = np.linalg.inv(M)
C = det_M * inv_M.T

# Compute antisymmetric part of C
K_cof = 0.5 * (C - C.T)

# --- Step 4: Calculate the final norm ---

# T is tridiagonal and orthogonally similar to K_cof.
# T^2 is orthogonally similar to K_cof^2. They have the same eigenvalues.
# K_cof is skew-symmetric, so K_cof^2 is symmetric.
K_cof_sq = K_cof @ K_cof

# For a symmetric matrix, singular values are the absolute values of its eigenvalues.
eigvals_K_cof_sq = np.linalg.eigvalsh(K_cof_sq)
singular_values = np.abs(eigvals_K_cof_sq)

# Sort singular values in descending order
singular_values = np.sort(singular_values)[::-1]

# The largest Ky Fan norm is the nuclear norm (k=rank), which is the sum of all singular values.
largest_ky_fan_norm = np.sum(singular_values)

# Final output: print the summation of singular values
# to show how the final result is obtained.
equation_parts = [f"{val:.6f}" for val in singular_values]
equation_str = " + ".join(equation_parts)
print(f"The singular values of the squared tridiagonal matrix are: {', '.join(equation_parts)}")
print("The largest Ky Fan norm is the sum of these singular values:")
print(f"{equation_str} = {largest_ky_fan_norm:.6f}")

<<<3.0>>>