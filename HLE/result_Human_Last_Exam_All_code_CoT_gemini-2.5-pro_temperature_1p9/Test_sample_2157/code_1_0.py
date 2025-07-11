import numpy as np
from numpy.linalg import det, inv, svd

def get_mandelbrot_matrix(n):
    """
    Generates the Mandelbrot matrix M_n based on the Griffin & Kraaijevanger definition.
    """
    N = 2**(n + 1) - 1
    if N <= 0:
        return np.array([])
    
    M = np.zeros((N, N))
    
    # First column based on the vector v
    for k in range(1, n + 2):
        idx = 2**k - 1
        if idx <= N:
            M[idx - 1, 0] = 1
            
    # Superdiagonal is 1 (backward shift operator)
    for i in range(N - 1):
        M[i, i + 1] = 1
            
    return M

def find_n0():
    """
    Finds the value of n that minimizes the given function F(n).
    """
    f_values = {}
    print("Step 1 & 2: Finding n_0")
    for n in range(4): # Check for n = 0, 1, 2, 3
        N = 2**(n + 1) - 1
        if N <= 0:
            continue
        
        M_n = get_mandelbrot_matrix(n)
        S_n = 0.5 * (M_n + M_n.T)
        
        # Check for existence of LDL' by checking principal minors
        decom_exists = True
        for k in range(1, N):
            minor_det = det(S_n[:k, :k])
            if np.isclose(minor_det, 0):
                decom_exists = False
                break
        
        if not decom_exists:
            print(f"For n = {n}, LDL' decomposition is undefined. F({n}) is not calculated.")
            f_values[n] = float('inf')
            continue

        # LDL' exists, proceed to calculate F(n)
        D_n_diag = []
        prev_minor_det = 1.0
        for k in range(1, N + 1):
            minor_det = det(S_n[:k, :k])
            D_n_diag.append(minor_det / prev_minor_det)
            prev_minor_det = minor_det
        D_n_diag = np.array(D_n_diag)
        
        det_D_n = np.prod(D_n_diag)
        trace_D_n = np.sum(D_n_diag)
        
        # Take abs for robustness with fractional powers of negative numbers
        geom_mean_part = np.power(np.abs(det_D_n), 1.0 / N)
            
        f_n = trace_D_n * geom_mean_part
        f_values[n] = f_n
        print(f"For n = {n}, F({n}) = {f_n:.4f}")

    # The analysis shows F(0)=1, F(1)=0, and F(n) is undefined for n>=2.
    # Thus, the minimum is at n=1.
    n_0 = min(f_values, key=f_values.get)
    print(f"\nThe function F(n) is minimized at n_0 = {n_0}\n")
    return n_0

def main():
    """
    Main function to execute the plan.
    """
    # Steps 1 & 2: Determine n_0
    n_0 = find_n0()

    # Step 3: Perform calculations for n_0
    print(f"Step 3: Calculations for n_0 = {n_0}")
    M_n0 = get_mandelbrot_matrix(n_0)

    # Cofactor matrix C = (det(M) * M^-1)^T
    det_M_n0 = det(M_n0)
    Cofactor_M = (det_M_n0 * inv(M_n0)).T
    
    # Antisymmetric part A = 0.5 * (C - C^T)
    A = 0.5 * (Cofactor_M - Cofactor_M.T)

    # Tridiagonal part T
    T = np.diag(np.diag(A, -1), -1) + np.diag(np.diag(A, 0), 0) + np.diag(np.diag(A, 1), 1)

    # Square of T, let's call it B
    B = T @ T
    print(f"The matrix B = T^2 is:\n{np.round(B, 4)}\n")
    
    # Step 4: Find the largest Ky Fan norm
    print("Step 4: Find the largest Ky Fan norm of B")
    singular_values = svd(B, compute_uv=False)
    # Sort in descending order
    singular_values = np.sort(singular_values)[::-1]
    
    ky_fan_norms = [np.sum(singular_values[:k+1]) for k in range(len(singular_values))]
    
    largest_ky_fan_norm = np.max(ky_fan_norms)

    print("The singular values of B are:", ", ".join([f"{s:.4f}" for s in singular_values]))
    print("The Ky Fan k-norms are:")
    for k, norm_val in enumerate(ky_fan_norms):
        print(f"  k={k+1}: ||B||_({k+1}) = {norm_val:.4f}")
        
    print(f"\nThe largest Ky Fan norm is the nuclear norm, with value {largest_ky_fan_norm:.4f}")

    # Output the final equation as requested
    s_vals_str = " + ".join([f"{s:.5f}" for s in singular_values])
    final_eq_str = f"Final Equation: {s_vals_str} = {largest_ky_fan_norm:.5f}"
    print("\n" + final_eq_str)
    
    return largest_ky_fan_norm

if __name__ == '__main__':
    result = main()
    # To conform to the specified output format
    # print(f"\n<<<{result:.5f}>>>")

main()