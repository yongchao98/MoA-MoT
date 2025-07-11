import numpy as np
from scipy.linalg import ldl, hessenberg, svd

def construct_mandelbrot_matrix(n):
    """
    Constructs the Mandelbrot matrix M_n of size (2^(n+1)-1) x (2^(n+1)-1).
    """
    N = 2**(n + 1) - 1
    M = np.zeros((N, N))
    
    if N > 1:
        # Set superdiagonal to 1
        np.fill_diagonal(M[:, 1:], 1)

    # Set the first column based on the standard definition
    for k in range(1, N + 1):
        i = k - 1
        if k == 1:
            M[i, 0] = -1
            continue
        
        # Binary representation of k, reversed for b_0, b_1, ... order
        b = bin(k)[2:][::-1]
        prod = 1
        for j in range(len(b) - 1):
            b_j = int(b[j])
            b_j_plus_1 = int(b[j+1])
            prod *= (-1)**(b_j * b_j_plus_1)
        M[i, 0] = -prod
        
    return M

def solve():
    """
    Main function to solve the entire problem.
    """
    # Part 1: Find n_0
    results = []
    print("Step 1: Finding n_0 by minimizing Tr(D_n) * (Det(D_n))^(1/n)")
    for n_test in range(1, 9):
        try:
            M_n = construct_mandelbrot_matrix(n_test)
            S_n = (M_n + M_n.T) / 2
            
            # Use LDL' decomposition for indefinite matrices
            # P S_n P.T = L D L.T, where D is block-diagonal
            lu, d, perm = ldl(S_n)
            
            trace_D = np.trace(d)
            det_D = np.linalg.det(d)

            # Avoid complex numbers from even roots of negative numbers
            if det_D < 0 and n_test % 2 == 0:
                print(f"  n={n_test}: Skipping (Det(D_n) < 0 for even n)")
                continue
            
            # Calculate value using real nth root
            term = np.sign(det_D) * (np.abs(det_D)**(1/n_test))
            val = trace_D * term
            
            results.append({'n': n_test, 'value': val})
            print(f"  n={n_test}, Value = {val:.4f}")

        except np.linalg.LinAlgError:
            print(f"  n={n_test}: Skipping (Singular matrix encountered)")
    
    if not results:
        print("Could not find a suitable n_0. Aborting.")
        return

    min_result = min(results, key=lambda x: x['value'])
    n0 = min_result['n']
    print(f"\nMinimum value found at n = {n0}. Thus, n_0 = {n0}.\n")

    # Part 2: Matrix operations for n_0
    print(f"Step 2: Performing calculations for n_0 = {n0}")
    M = construct_mandelbrot_matrix(n0)
    N = M.shape[0]
    print(f"  Constructed M_{n0} of size {N}x{N}.")

    # Cofactor matrix calculation
    det_M = np.linalg.det(M)
    inv_M = np.linalg.inv(M)
    cof_M = det_M * inv_M.T
    print(f"  Calculated cofactor matrix of M_{n0}. Det(M_{n0}) = {det_M:.2f}.")
    
    # Antisymmetric part
    A = (cof_M - cof_M.T) / 2
    print("  Calculated the antisymmetric part (A) of the cofactor matrix.")

    # Tridiagonalization
    T_A = hessenberg(A)
    # For a skew-symmetric matrix, the Hessenberg form is skew-tridiagonal.
    print("  Reduced A to a skew-tridiagonal matrix (T_A).")
    
    # Square the matrix
    T_A_sq = T_A @ T_A
    print("  Calculated the square of the tridiagonal matrix (T_A_sq).")

    # Part 3: Compute the largest Ky Fan norm (Nuclear Norm)
    print("\nStep 3: Calculating the largest Ky Fan norm (Nuclear Norm) of T_A_sq")
    s = svd(T_A_sq, compute_uv=False)
    nuclear_norm = np.sum(s)

    # Output the final equation as requested
    s_str = " + ".join([f"{val:.4f}" for val in s])
    
    print(f"\nThe singular values of T_A_sq are:\n{s}")
    print("\nThe final equation for the nuclear norm is:")
    # For readability, print in chunks if too long
    if len(s_str) > 120:
        s_chunks = [s_str[i:i+120] for i in range(0, len(s_str), 120)]
        print(" + \n".join(s_chunks))
    else:
        print(s_str)

    print(f"\nResult of the sum: {nuclear_norm:.4f}")
    
    print("\nThe largest Ky Fan norm is the sum of these singular values.")
    print(f"Final Answer: {nuclear_norm}")
    print(f"\n<<<{nuclear_norm:.8f}>>>")

if __name__ == '__main__':
    solve()