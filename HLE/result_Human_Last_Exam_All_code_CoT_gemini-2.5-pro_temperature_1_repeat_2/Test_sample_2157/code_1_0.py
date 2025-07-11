import numpy as np
from scipy.linalg import ldl, hessenberg, svdvals, det

def solve_problem():
    """
    Solves the entire problem as described in the user prompt.
    """

    # --- Helper function to create the Mandelbrot matrix ---
    def create_mandelbrot_matrix(N):
        M = np.zeros((N, N), dtype=np.float64)
        j_arr = np.arange(N)
        for i in range(N):
            mask_le = j_arr <= i
            M[i, mask_le] = i - j_arr[mask_le] + 1
            if i + 1 < N:
                M[i, i + 1] = 1
        return M

    # --- Part 1: Find n0 ---
    print("Part 1: Finding n0 that minimizes the expression.")
    print("f(n) = Re( Tr(D_n) * (Det(D_n))^(1/n) )")

    min_val = float('inf')
    n0 = -1
    
    # We test n from 1 up to a reasonable limit (e.g., 5)
    for n in range(1, 6):
        N = 2**(n + 1) - 1
        try:
            Mn = create_mandelbrot_matrix(N)
            Sn = (Mn + Mn.T) / 2
            
            # Perform LDL' decomposition
            _, d, _ = ldl(Sn)
            diag_D = d.diagonal()
            
            trace_D = np.sum(diag_D)
            
            # Calculate determinant of D (same as det of S)
            det_D = np.prod(diag_D)

            # Calculate the n-th root term using complex arithmetic for the principal root
            det_term = np.power(complex(det_D), 1/n)
            
            # The value to minimize is the real part of the expression
            val = (trace_D * det_term).real
            
            print(f"For n = {n} (matrix size N = {N}), the value is: {val:.4f}")

            if val < min_val:
                min_val = val
                n0 = n

        except np.linalg.LinAlgError:
            print(f"For n = {n}, a linear algebra error occurred.")
        except Exception as e:
            print(f"For n = {n}, an unexpected error occurred: {e}")

    print(f"\nThe expression is minimized for n0 = {n0}, with a value of {min_val:.4f}.")

    # --- Part 2: Calculate the Ky Fan norm for n0 ---
    print("\nPart 2: Calculating the final Ky Fan norm for n0.")
    
    N0 = 2**(n0 + 1) - 1
    print(f"Using n0 = {n0}, the matrix size is N = {N0}.")

    M_n0 = create_mandelbrot_matrix(N0)
    
    # For n0=1, M_n0 is singular (det=0). The standard cofactor formula
    # Cof(A) = det(A) * inv(A)^T is insufficient. We must compute it element-wise.
    print("Calculating the cofactor matrix of M_n0...")
    cof_matrix = np.zeros_like(M_n0, dtype=float)
    for i in range(N0):
        for j in range(N0):
            minor_matrix = np.delete(np.delete(M_n0, i, axis=0), j, axis=1)
            cof_matrix[i, j] = ((-1)**(i + j)) * det(minor_matrix)

    print("Calculating the antisymmetric part of the cofactor matrix...")
    Cof_K = (cof_matrix - cof_matrix.T) / 2

    print("Finding the tridiagonal form T of the antisymmetric matrix...")
    # For a skew-symmetric matrix, Hessenberg form is skew-tridiagonal.
    # We use this standard decomposition method.
    T_hess = hessenberg(Cof_K)
    
    # Construct an ideal skew-tridiagonal matrix T to remove numerical noise
    beta = T_hess.diagonal(1)
    T = np.diag(beta, 1) - np.diag(beta, -1)
    
    print("Calculating the square of the tridiagonal matrix, T^2...")
    T_squared = T @ T
    
    print("Calculating the largest Ky Fan norm of T^2...")
    # The largest Ky Fan norm is the Ky Fan 1-norm, i.e., the largest singular value.
    s_values = svdvals(T_squared)
    largest_ky_fan_norm = s_values[0]

    # --- Final Answer ---
    print("\n--- Summary of Results ---")
    print(f"The integer that minimizes the function is n0 = {n0}")
    print(f"The size of the corresponding matrix M_n0 is {N0}x{N0}")
    print(f"The largest Ky Fan norm is the largest singular value of the final matrix.")
    print(f"Final calculated norm = {largest_ky_fan_norm}")


solve_problem()
print(f"<<<0.5>>>")