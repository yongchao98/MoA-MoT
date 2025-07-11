import numpy as np

def build_mandelbrot_matrix(n):
    """
    Constructs the Mandelbrot matrix M_n of size (2^(n+1)-1) x (2^(n+1)-1).
    """
    N = 2**(n + 1) - 1
    M = np.zeros((N, N))
    
    # Diagonal and superdiagonal
    for i in range(N):
        M[i, i] = -1
        if i < N - 1:
            M[i, i + 1] = 1
            
    # Subdiagonal
    subdiagonal_zeros_at_i = {2**k - 1 for k in range(2, n + 2)}
    for i in range(N - 1):
        if (i + 1) not in subdiagonal_zeros_at_i:
            M[i + 1, i] = -1
            
    return M

def get_ldl_diag(S):
    """
    Performs LDL' decomposition and returns the diagonal matrix D.
    This is a simple implementation for the specific block-diagonal structure.
    """
    N = S.shape[0]
    D_diag = np.zeros(N)
    i = 0
    while i < N:
        if i + 1 < N and S[i, i + 1] != 0:
            # 2x2 block
            d1 = S[i, i]
            l21 = S[i+1, i] / d1
            d2 = S[i+1, i+1] - l21**2 * d1
            D_diag[i] = d1
            D_diag[i+1] = d2
            i += 2
        else:
            # 1x1 block
            D_diag[i] = S[i, i]
            i += 1
    return D_diag

def calculate_minimization_function(n):
    """
    Calculates the value of the function to be minimized for a given n.
    """
    Mn = build_mandelbrot_matrix(n)
    Sn = 0.5 * (Mn + Mn.T)
    
    # We use -Sn for LDL' decomposition as it's positive definite
    neg_Sn = -Sn
    
    # D_n is the diagonal of the LDL' decomposition
    D_n_diag = get_ldl_diag(neg_Sn)
    
    trace_Dn = np.sum(D_n_diag)
    det_Dn = np.prod(D_n_diag)
    
    # Calculate f(n)
    f_n = trace_Dn * (det_Dn)**(1/n)
    return f_n

def solve():
    """
    Main function to solve the problem.
    """
    # Step 1 & 2: Find n_0 by observing the trend of the minimization function
    print("Step 1: Find n_0 by evaluating the function f(n) for small n.")
    f_values = []
    for n in range(1, 5):
        f_n = calculate_minimization_function(n)
        f_values.append(f_n)
        print(f"f({n}) = {f_n:.4f}")
    
    # The function is increasing, so the minimum is at n_0=1
    n0 = np.argmin(f_values) + 1
    print(f"\nThe function is minimized at n_0 = {n0}.\n")

    # Step 3: Matrix calculations for n_0 = 1
    print(f"Step 2: Perform matrix calculations for n_0 = {n0}.")
    M1 = build_mandelbrot_matrix(n0)
    print("M_1:\n", M1)

    # Calculate cofactor matrix C of M1
    # C = (inv(M1).T) * det(M1)
    det_M1 = np.linalg.det(M1)
    C = np.linalg.inv(M1).T * det_M1
    print("\nCofactor matrix C:\n", np.round(C))

    # Calculate the antisymmetric part A
    A = 0.5 * (C - C.T)
    print("\nAntisymmetric part A:\n", np.round(A))

    # Step 4: Parlett-Reid, T^2, and Ky Fan Norm
    # Since A is already tridiagonal, the tridiagonal matrix T is A itself.
    T = A
    print("\nThe matrix A is already tridiagonal, so T = A.")

    # Calculate the square of T
    B = T @ T
    print("\nB = T^2:\n", np.round(B))

    # Singular values of B are the absolute values of its eigenvalues (since B is symmetric)
    eigenvalues = np.linalg.eigvals(B)
    singular_values = np.sort(np.abs(np.round(eigenvalues, 5)))[::-1]
    print(f"\nThe singular values of B are: {singular_values[0]}, {singular_values[1]}, {singular_values[2]}")

    # The largest Ky Fan norm is the nuclear norm (k=N), sum of all singular values.
    ky_fan_norms = [np.sum(singular_values[:k+1]) for k in range(len(singular_values))]
    largest_norm = np.max(ky_fan_norms)
    
    print("\nThe Ky Fan k-norms are:")
    print(f"k=1: {ky_fan_norms[0]}")
    print(f"k=2: {ky_fan_norms[1]}")
    print(f"k=3: {ky_fan_norms[2]}")

    print(f"\nThe largest Ky Fan norm is the sum of all singular values:")
    print(f"{singular_values[0]} + {singular_values[1]} + {singular_values[2]} = {largest_norm}")
    
    print("\nFinal Answer:")
    print(f"<<<{largest_norm}>>>")

solve()