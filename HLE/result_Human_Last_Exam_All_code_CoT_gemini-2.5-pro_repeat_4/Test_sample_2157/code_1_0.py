import numpy as np
from scipy.linalg import ldl, inv, hessenberg, svdvals

def construct_M(n):
    """Constructs the Mandelbrot Matrix M_n of order (2**(n+1)-1)."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    N = 2**(n + 1) - 1
    M = np.zeros((N, N))
    
    # Using 0-based indexing for the matrix
    # Rule 1: M[i, i+1] = 1
    for i in range(N - 1):
        M[i, i + 1] = 1
        
    # Rule 2: M[i-1, 0] = -1 if i is even (for 1-based index i)
    for i_one_based in range(2, N + 1):
        if i_one_based % 2 == 0:
            M[i_one_based - 1, 0] = -1
            
    # Rule 3: M[i-1, j-1] = -1 if i = 2j - 1 for j >= 2 (1-based indices)
    for j_one_based in range(2, N + 1):
        i_one_based = 2 * j_one_based - 1
        if i_one_based <= N:
            M[i_one_based - 1, j_one_based - 1] = -1
            
    return M

def calculate_objective_function(n, M):
    """Calculates f(n) = Tr(D_n) * (Det(D_n))**(1/n)."""
    # Assumption: "symmetric part" refers to the positive semi-definite matrix M @ M.T
    S = M @ M.T
    
    # The determinant of M_n is 0 for odd n and 1 for even n.
    # So, det(S) = det(M)^2 is 0 for odd n and 1 for even n.
    # This means det(D) from LDL' will also be 0 for odd n.
    det_M = np.linalg.det(M)
    if abs(det_M) < 1e-9:
        return 0.0

    # For even n, S is positive definite, so LDL' is stable and D is diagonal.
    try:
        _, d, _ = ldl(S)
        D_diag = d.diagonal()
        trace_D = np.sum(D_diag)
        det_D = np.prod(D_diag)
        
        # For even n, det(D) should be close to 1.
        # For n=2, det(M) is not 1, so we calculate it.
        val = trace_D * (np.abs(det_D))**(1/n)
        return val
    except np.linalg.LinAlgError:
        return None

def main():
    """Main function to solve the problem."""
    print("Step 1: Find the value of n_0 that minimizes the objective function.")
    print("Using the interpretation that the 'symmetric part' of M_n is M_n @ M_n.T\n")

    results = {}
    for n in range(1, 5):
        M = construct_M(n)
        f_n = calculate_objective_function(n, M)
        results[n] = f_n
        print(f"For n = {n}, f(n) = {f_n}")

    min_val = min(results.values())
    minimizers = [n for n, val in results.items() if val == min_val]
    
    print(f"\nThe minimum value of the function is {min_val}, achieved for n = {minimizers}.")
    
    # Choose the simplest minimizer
    n_0 = minimizers[0]
    print(f"We choose the simplest minimizer, n_0 = {n_0}.\n")
    
    print(f"Step 2: Calculate the Ky Fan norm for n_0 = {n_0}.")
    
    # For n_0 = 1, the calculation is straightforward.
    if n_0 == 1:
        M1 = construct_M(1)
        print("M_1 matrix:\n", M1)
        
        # The cofactor matrix C_1 of M_1 can be computed analytically.
        # C_ij = (-1)**(i+j) * det(M_ij)
        C1 = np.array([[1., 0., 1.],
                       [0., 0., 0.],
                       [1., 0., 1.]])
        print("\nCofactor matrix C_1:\n", C1)
        
        # The antisymmetric part A_1 = 0.5 * (C_1 - C_1.T)
        A1 = 0.5 * (C1 - C1.T)
        print("\nAntisymmetric part A_1:\n", A1)
        
        # The tridiagonalization of a zero matrix is the zero matrix.
        T = np.zeros_like(A1)
        print("\nTridiagonal matrix T of A_1:\n", T)
        
        # Square of the tridiagonal matrix
        T2 = T @ T
        print("\nSquare of the tridiagonal matrix, T^2:\n", T2)
        
        # The Ky Fan 1-norm is the largest singular value.
        # The singular values of a zero matrix are all zero.
        final_norm = 0.0
        
        print("\nThe final equation is: max(svd(T^2))")
        print(f"T^2 is the zero matrix, so its singular values are all 0.")
        print(f"The largest singular value is {final_norm}.")

    else:
        # This case is not expected based on the analysis.
        print("Analysis for n_0 > 1 is more complex and not required as n_0=1 is a minimizer.")

if __name__ == '__main__':
    main()
