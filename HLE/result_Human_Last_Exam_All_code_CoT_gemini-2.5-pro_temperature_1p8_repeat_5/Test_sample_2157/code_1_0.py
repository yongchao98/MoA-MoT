import numpy as np
from numpy.linalg import inv, det

def solve_mandelbrot_matrix_problem():
    """
    Solves the complex matrix problem described by the user.
    The solution is based on interpreting the Mandelbrot matrix M_n
    as the companion matrix for a specific polynomial related to the
    Mandelbrot set definition. We proceed assuming n_0 = 1, which corresponds
    to a 3x3 matrix.
    """

    # Step 1: Define the matrix M_1 for n=1
    # Polynomial p_2(z)/z = z^3 + 2z^2 + z + 1
    # Companion matrix M_1 (in upper Hessenberg form)
    M1 = np.array([
        [0., 1., 0.],
        [0., 0., 1.],
        [-1., -1., -2.]
    ])
    print(f"Using M_1 (for n=1):\n{M1}\n")

    # Step 2: Compute the cofactor matrix C_1
    # C = det(M) * (M^T)^-1
    det_M1 = det(M1)
    # Note: A small tolerance is used because numerical precision might make det not exactly -1.
    if not np.isclose(det_M1, -1.0):
        print(f"Warning: Determinant is {det_M1}, not exactly -1. Using precise value.")
    
    C1 = det_M1 * inv(M1.T)
    print(f"Determinant of M_1: {det_M1:.4f}")
    print(f"Cofactor matrix C_1:\n{C1}\n")

    # Step 3: Compute the antisymmetric part K'_1
    K_prime_1 = 0.5 * (C1 - C1.T)
    print(f"Antisymmetric part of the cofactor matrix, K'_1:\n{K_prime_1}\n")

    # Step 4: Find the eigenvalues of K'_1
    # For a real skew-symmetric matrix, eigenvalues are purely imaginary.
    eigenvalues = np.linalg.eigvals(K_prime_1)
    print(f"Eigenvalues of K'_1: [{', '.join(f'{v:.4f}' for v in eigenvalues)}]")
    
    # The tridiagonal matrix T has the same eigenvalues.
    # The spectral norm ||K'_1||_2 is the max absolute value of its eigenvalues.
    spectral_norm_K_prime = np.max(np.abs(eigenvalues))
    print(f"The largest Ky Fan norm of T (spectral norm) is max|lambda| = {spectral_norm_K_prime:.4f}")

    # Step 5: Compute the largest Ky Fan norm of T^2
    # This is equivalent to (||K'_1||_2)^2
    result = spectral_norm_K_prime**2
    
    # Final equation and result
    print(f"\nThe final value is the square of this norm.")
    print(f"Equation: ({spectral_norm_K_prime:.4f})^2 = {result:.4f}")
    
    # Print the exact fraction for clarity
    # The characteristic polynomial of K'_1 is lambda^3 + 2.75 * lambda = 0
    # The non-zero eigenvalues are +/- i*sqrt(2.75)
    # The squared norm is 2.75, or 11/4
    print(f"The exact value is 11/4.")
    print(f"Final answer: {result}")
    

solve_mandelbrot_matrix_problem()
<<<2.75>>>