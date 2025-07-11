import numpy as np
from scipy.linalg import ldl
import cmath

def get_mandelbrot_matrix(n):
    """
    Constructs the matrix M_n of size (2n+1) x (2n+1).
    This family is sparse, upper Hessenberg, and matches the problem's size
    formula for n=1.
    """
    if n == 0:
        return np.array([[-1.0]])
    
    m_prev = get_mandelbrot_matrix(n - 1)
    prev_size = m_prev.shape[0]
    size = prev_size + 2
    
    m = np.zeros((size, size))
    m[:prev_size, :prev_size] = m_prev
    
    m[0, prev_size] = 1.0
    m[0, prev_size + 1] = 1.0
    
    m[prev_size, prev_size - 1] = 1.0
    m[prev_size, prev_size + 1] = -1.0
    
    m[prev_size + 1, prev_size] = 1.0
    m[prev_size + 1, prev_size + 1] = 1.0
    
    return m

def get_ldl_diagonal(matrix):
    """
    Computes the diagonal matrix D of the LDL' decomposition.
    Returns the diagonal elements as a numpy array.
    This works if all leading principal minors are non-zero.
    """
    n = matrix.shape[0]
    d = np.zeros(n)
    # Using a temporary matrix to avoid modifying the original
    a = np.copy(matrix)
    for i in range(n):
        # Check if pivot is zero
        if np.abs(a[i, i]) < 1e-12:
            return None # Decomposition fails without pivoting
        d[i] = a[i, i]
        for j in range(i + 1, n):
            factor = a[j, i] / d[i]
            a[j, i:] -= factor * a[i, i:]
    return d

def solve():
    """
    Main function to solve the entire problem.
    """
    print("Step 1: Find the value of n_0 that minimizes F(n).")
    f_values = []
    min_val = float('inf')
    n0 = -1

    for n in range(1, 6):
        print(f"\n--- Calculating for n = {n} ---")
        M = get_mandelbrot_matrix(n)
        S = 0.5 * (M + M.T)
        
        # Check if LDL' is possible
        minors_ok = True
        for i in range(1, S.shape[0] + 1):
            if np.abs(np.linalg.det(S[:i, :i])) < 1e-9:
                minors_ok = False
                break
        
        if not minors_ok:
            print(f"LDL' decomposition without pivoting is not possible for n={n}.")
            continue

        D_diag = get_ldl_diagonal(S)
        
        trace = np.sum(D_diag)
        det = np.prod(D_diag)
        
        # Use complex power for robustness with negative determinants
        f_n = trace * cmath.exp((1/n) * cmath.log(det))
        
        # We are minimizing a real value, so we take the real part.
        f_n_real = f_n.real
        
        f_values.append((n, f_n_real))
        print(f"M_n is a {M.shape[0]}x{M.shape[0]} matrix.")
        print(f"Diagonal of D_n: {np.round(D_diag, 4)}")
        print(f"Tr(D_n) = {trace:.4f}")
        print(f"Det(D_n) = {det:.4f}")
        print(f"F(n) = Tr(D_n) * (Det(D_n))^(1/n) = {f_n_real:.4f}")

        if f_n_real < min_val:
            min_val = f_n_real
            n0 = n

    print("\n-----------------------------------------")
    print(f"F(n) values: {[(n, round(v, 4)) for n, v in f_values]}")
    print(f"The minimum value of F(n) occurs at n_0 = {n0}")
    print("-----------------------------------------")

    print(f"\nStep 2: Calculate the final Ky Fan norm for n_0 = {n0}.")
    M_n0 = get_mandelbrot_matrix(n0)
    print(f"M_{n0} is a {M_n0.shape[0]}x{M_n0.shape[0]} matrix.")

    det_M_n0 = np.linalg.det(M_n0)
    if abs(det_M_n0) < 1e-9:
        print("Determinant is zero, cannot proceed.")
        return

    inv_M_n0 = np.linalg.inv(M_n0)
    
    # Antisymmetric part of the cofactor matrix
    # A = det(M)/2 * ((M^-1)^T - M^-1)
    # We need the norm of A^2, which is ||A||^2 for a normal matrix A.
    
    K_skew = 0.5 * (inv_M_n0.T - inv_M_n0)
    
    # Spectral norm is the largest singular value
    norm_K_skew = np.linalg.norm(K_skew, 2)
    
    norm_A_n0_sq = (det_M_n0**2) * (norm_K_skew**2)
    
    print(f"Determinant of M_{n0}: det = {det_M_n0:.4f}")
    print(f"Inverse of M_{n0} calculated.")
    print(f"Skew-symmetric part of the inverse is K = (M_{n0}^-1 - (M_{n0}^-1)^T)/2.")
    print(f"Spectral norm of K: ||K||_2 = {norm_K_skew:.4f}")
    print(f"The matrix A_{n0} is det * K.")
    print(f"The required value is the largest Ky Fan norm of T_{n0}^2, which is ||A_{n0}||_2^2.")
    print(f"Final calculation: ||A_{n0}||_2^2 = (det)^2 * ||K||_2^2")
    print(f"Result = ({det_M_n0:.4f})^2 * ({norm_K_skew:.4f})^2 = {norm_A_n0_sq:.4f}")
    
    print(f"\nFinal Answer: {norm_A_n0_sq}")
    
solve()
<<<1.5>>>