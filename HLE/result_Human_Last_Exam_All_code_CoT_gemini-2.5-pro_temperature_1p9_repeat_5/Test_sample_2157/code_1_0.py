import numpy as np

def solve_mandelbrot_matrix_problem():
    """
    Solves the complex matrix problem described by the user.
    The solution follows a multi-step derivation:
    1. A specific family of "Mandelbrot Matrices" Mn is chosen that satisfies
       all given properties (sparse, upper Hessenberg, eigenvalues on the
       Mandelbrot set boundary). The chosen Mn is a real, block-diagonal matrix
       with eigenvalues {i, -i, -2}, e.g., for n=1, a 3x3 matrix.
    2. For this Mn, its symmetric part Sn is singular. This means the
       minimization criterion Tr(Dn) * (Det(Dn))^(1/n) equals 0 for any n >= 1.
       This makes the choice of n_0 arbitrary among values >= 1.
       We choose n_0 = 1 for demonstration. The final result is independent of this choice.
    3. The rest of the problem is a direct calculation using n_0=1.
    """
    
    n0 = 1
    N = 2**(n0 + 1) - 1 # Matrix size for n=1 is 3x3

    # Step 1: Define the Matrix M_{n_0} for n_0 = 1
    # This matrix is tridiagonal (a form of sparse upper Hessenberg)
    # and has eigenvalues {i, -i, -2}, which lie on the Mandelbrot set boundary.
    M_n0 = np.array([
        [0., -1., 0.],
        [1.,  0., 0.],
        [0.,  0., -2.]
    ])
    
    print(f"Step 1: The chosen Mandelbrot Matrix M_{n0} for n0={n0} (size {N}x{N}) is:")
    print(M_n0)
    print("-" * 20)
    
    # Step 2: Compute the cofactor matrix of M_{n_0}
    # For a 3x3 matrix A = [[a,b,c],[d,e,f],[g,h,i]], the cofactor matrix C is
    # [[(ei-fh), -(di-fg), (dh-eg)],
    #  [-(bi-ch), (ai-cg), -(ah-bg)],
    #  [(bf-ce), -(af-cd), (ae-bd)]]
    # Alternatively, C = (det(A) * inv(A)^T)
    det_M = np.linalg.det(M_n0)
    # Adding a small epsilon to diagonal for stability if matrix is singular
    # but in our case det is -2 so it's not singular.
    inv_M_T = np.linalg.inv(M_n0).T
    C_n0 = det_M * inv_M_T

    print(f"Step 2: The cofactor matrix C_{n0} of M_{n0} is:")
    print(C_n0)
    print("-" * 20)
    
    # Step 3: Compute the antisymmetric part of the cofactor matrix
    A_prime_n0 = (C_n0 - C_n0.T) / 2
    
    print(f"Step 3: The antisymmetric part of the cofactor matrix, A'_{n0}, is:")
    print(A_prime_n0)
    print("-" * 20)

    # Step 4: Parlett-Reid decomposition of A'_{n_0} to get a tridiagonal matrix T_{n_0}
    # The Parlett-Reid algorithm finds a permutation P such that P*A*P^T is tridiagonal.
    # In our case, A'_{n_0} is already tridiagonal, so the permutation is the identity matrix,
    # and the tridiagonal matrix T_{n_0} is A'_{n_0} itself.
    T_n0 = A_prime_n0
    
    print("Step 4: The tridiagonal matrix T_{n0} from Parlett-Reid decomposition is:")
    print(T_n0)
    print("-" * 20)

    # Step 5: Compute the square of T_{n_0}
    T_n0_squared = T_n0 @ T_n0
    
    print(f"Step 5: The square of the tridiagonal matrix, T_{n0}^2, is:")
    print(T_n0_squared)
    print("-" * 20)

    # Step 6: Find the largest Ky Fan norm of T_{n_0}^2
    # The largest Ky Fan norm is the Ky Fan 1-norm, which is the spectral norm.
    # The spectral norm is the largest singular value of the matrix.
    # For a normal matrix like T_n0_squared, singular values are the absolute values of eigenvalues.
    singular_values = np.linalg.svd(T_n0_squared, compute_uv=False)
    largest_ky_fan_norm = singular_values[0]

    # "output each number in the final equation!"
    # The final equation is ||T_{n_0}^2||_2 = result
    print("Step 6: The final calculation is the largest Ky Fan norm (spectral norm).")
    print(f"The singular values of T_{n0}^2 are: {singular_values}")
    print(f"The largest singular value is: {largest_ky_fan_norm:.2f}")
    
    # Final Result
    final_answer = np.round(largest_ky_fan_norm)
    print("\nFinal Answer:")
    print(final_answer)
    

solve_mandelbrot_matrix_problem()