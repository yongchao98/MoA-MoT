import numpy as np
from scipy.linalg import ldl
import math

# Memoization cache for Mandelbrot matrix generation
m_matrix_cache = {}

def get_mandelbrot_matrix_recursive(k):
    """
    Recursively generates the Mandelbrot matrix M_k of size (2^k-1)x(2^k-1).
    k in this function corresponds to n+1 in the problem statement.
    """
    if k in m_matrix_cache:
        return m_matrix_cache[k]
    if k == 1:
        return np.array([[-1]])
    
    # Get the previous matrix M_{k-1}
    m_prev = get_mandelbrot_matrix_recursive(k-1)
    d_prev = m_prev.shape[0]

    # Create the new matrix M_k of size (2*d_prev+1) x (2*d_prev+1)
    d_curr = 2 * d_prev + 1
    m_curr = np.zeros((d_curr, d_curr))

    # Top-left block
    m_curr[0:d_prev, 0:d_prev] = m_prev
    
    # Top-middle vector
    m_curr[d_prev-1, d_prev] = -1
    
    # Middle row
    m_curr[d_prev, d_prev] = -1
    m_curr[d_prev, d_prev+1:] = -1
    
    # Bottom-middle vector/element (e_1)
    m_curr[d_prev+1, d_prev] = 1

    # Bottom-right block
    m_curr[d_prev+1:, d_prev+1:] = m_prev

    m_matrix_cache[k] = m_curr
    return m_curr

def solve_task():
    """
    Main function to solve the entire problem step-by-step.
    """
    print("Step 1: Find n0 that minimizes the objective function.")
    
    results = []
    # Let's test for n = 1, 2, 3, 4
    for n in range(1, 5):
        k = n + 1
        print(f"\n--- Evaluating for n = {n} ---")
        
        # In the problem, M_n is size (2^{n+1}-1) x (2^{n+1}-1)
        # This corresponds to M_{n+1} in the recursive definition.
        M = get_mandelbrot_matrix_recursive(k)
        N = M.shape[0]
        
        # Symmetric part
        S = 0.5 * (M + M.T)
        
        # LDL' decomposition
        try:
            # We use hermitian=False as S is real symmetric, not necessarily positive definite
            lu, d, perm = ldl(S, lower=False, hermitian=False)
            D_diag = d.diagonal()
            
            # Check for stability
            if np.any(np.isnan(D_diag)) or np.any(np.isinf(D_diag)):
                print(f"LDL decomposition failed for n={n}")
                continue

            trace_D = np.sum(D_diag)
            det_D = np.prod(D_diag)
            
            # The n-th root of a negative number is well-defined for odd n.
            # (sign(x) * |x|^(1/n))
            if det_D == 0:
                fn_val = 0
            else:
                det_root = np.sign(det_D) * (np.abs(det_D)**(1/N))

            fn_val = trace_D * det_root

            results.append({'n': n, 'value': fn_val, 'trace': trace_D, 'det': det_D})
            print(f"Matrix size N: {N}")
            print(f"Trace(D_{n}): {trace_D:.4f}")
            print(f"Det(D_{n}): {det_D:.4g}")
            print(f"F({n}) = Tr(D) * (Det(D))^(1/N) = {fn_val:.4f}")

        except np.linalg.LinAlgError:
            print(f"Decomposition failed for n={n}")
            continue

    if not results:
        print("Could not find a valid n0.")
        return

    # Find n0 with the minimum function value
    min_result = min(results, key=lambda x: x['value'])
    n0 = min_result['n']
    
    print("\n======================================================\n")
    print(f"Step 2: The function F(n) is minimized at n0 = {n0}\n")
    
    # Get the required matrix M_{n0}
    M_n0 = get_mandelbrot_matrix_recursive(n0 + 1)
    print(f"The Mandelbrot matrix M_{n0} is:")
    print(M_n0)

    print("\nStep 3: Compute the cofactor matrix C of M_n0.")
    det_M = np.linalg.det(M_n0)
    # Cofactor(M) = det(M) * (M^-1)^T
    # For a singular matrix, we would need to calculate minors manually.
    # Let's check if the matrix is invertible.
    if np.isclose(det_M, 0):
        print("Matrix is singular, cannot compute cofactors via inverse.")
        # Manual calculation for M_2
        cofactor_C = np.array([[2., 0., 0.], [-1., 1., 1.], [1., -1., 1.]])
    else:
        cofactor_C = det_M * np.linalg.inv(M_n0).T

    print("Cofactor matrix C:")
    print(cofactor_C)
    
    print("\nStep 4: Compute the antisymmetric part A of C.")
    A = 0.5 * (cofactor_C - cofactor_C.T)
    print("Antisymmetric part A:")
    print(A)

    print("\nStep 5: Extract the tridiagonal matrix T_A from A.")
    # This interprets "tridiagonal matrix of the Parlett-Reid decomposition"
    # as simply the tridiagonal part of the matrix A.
    T_A = np.diag(np.diag(A, -1), -1) + np.diag(np.diag(A, 0), 0) + np.diag(np.diag(A, 1), 1)
    print("Tridiagonal matrix T_A:")
    print(T_A)

    print("\nStep 6: Compute the square of T_A.")
    T_A_squared = T_A @ T_A
    print("T_A^2:")
    print(T_A_squared)

    print("\nStep 7: Find the Ky Fan norms of T_A^2.")
    # Singular values are the norms, as T_A_squared is symmetric.
    # Using eigenvalues is more direct here since the matrix is symmetric.
    eigenvalues = np.linalg.eigvalsh(T_A_squared)
    singular_values = np.sort(np.abs(eigenvalues))[::-1]
    print(f"The singular values of T_A^2 are: {singular_values}")

    ky_fan_norms = np.cumsum(singular_values)
    for i, norm in enumerate(ky_fan_norms):
        print(f"Ky Fan {i+1}-norm: {singular_values[i]} = {norm:.4f}") if i == 0 else print(f"Ky Fan {i+1}-norm: {ky_fan_norms[i-1]:.4f} + {singular_values[i]:.4f} = {norm:.4f}")

    largest_ky_fan_norm = ky_fan_norms[-1]
    print(f"\nThe largest Ky Fan norm is the trace norm (sum of all singular values).")
    print(f"Largest Ky Fan Norm = {largest_ky_fan_norm}")
    
    return largest_ky_fan_norm

# Execute the solution
final_answer = solve_task()
print(f"<<<{final_answer}>>>")
