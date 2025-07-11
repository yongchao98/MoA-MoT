import numpy as np

def solve():
    """
    Solves the matrix problem based on the step-by-step plan.
    """
    # Step 1 & 2: Based on analysis, we use n_0=1 and construct M_1.
    # The matrix M_1 is 3x3, sparse, upper Hessenberg, S_1=I,
    # and has eigenvalues {1, 1+i, 1-i} which lie on the Mandelbrot boundary.
    # The simplest such matrix is:
    M1 = np.array([
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0]
    ])

    # Step 3: Compute the cofactor matrix of M1.
    # For any invertible matrix A, adj(A) = det(A) * inv(A).
    # The cofactor matrix Cof(A) is the transpose of the adjugate matrix adj(A).
    # So, Cof(M1) = (det(M1) * inv(M1))^T
    det_M1 = np.linalg.det(M1)
    inv_M1 = np.linalg.inv(M1)
    adj_M1 = det_M1 * inv_M1
    Cof_M1 = adj_M1.T

    # Compute the antisymmetric part of the cofactor matrix.
    A1 = 0.5 * (Cof_M1 - Cof_M1.T)

    # A1 is already tridiagonal, so we set T1 = A1.
    T1 = A1
    
    # Step 4: Compute the square of T1.
    T1_squared = T1 @ T1
    
    # Compute the singular values of T1_squared.
    # Since T1_squared is a real symmetric matrix, its singular values
    # are the absolute values of its eigenvalues.
    eigenvalues = np.linalg.eigvalsh(T1_squared)
    
    # Sort singular values in descending order.
    singular_values = np.sort(np.abs(eigenvalues))[::-1]

    # The Ky Fan k-norm is the sum of the k largest singular values.
    # The "largest Ky Fan norm" is the nuclear norm, which is the sum
    # of all singular values (k=rank of the matrix).
    
    # Filter out near-zero singular values for clean printing
    significant_sv = [s for s in singular_values if not np.isclose(s, 0)]
    largest_ky_fan_norm = np.sum(significant_sv)
    
    print("Based on our analysis, we analyze the matrix M_1 for n_0=1:")
    print(M1)
    
    print("\nThe antisymmetric part of its cofactor matrix (T1) is:")
    print(T1)
    
    print("\nThe square of T1 is:")
    print(T1_squared)
    
    print("\nThe singular values of T1_squared are (sorted):")
    sv_str_list = [f"{s:.2f}" for s in singular_values]
    print(", ".join(sv_str_list))
    
    print("\nThe largest Ky Fan norm is the sum of the non-zero singular values.")
    sv_sum_str_list = [f"{s:.0f}" for s in significant_sv]
    final_eq = " + ".join(sv_sum_str_list)
    print(f"The equation for the largest Ky Fan norm is: {final_eq} = {largest_ky_fan_norm:.0f}")

solve()
<<<2>>>