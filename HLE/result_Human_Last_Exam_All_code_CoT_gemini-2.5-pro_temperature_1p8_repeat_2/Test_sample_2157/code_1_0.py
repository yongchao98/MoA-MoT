import numpy as np

def solve_matrix_problem():
    """
    This function solves the user's request by following the outlined plan.
    """
    # Step 1 & 2: Identify M_n and determine n_0.
    # We conclude n_0 = 1, which means the matrix size is N=3.
    # M_1 is a 3x3 Jordan block with eigenvalue -2.
    M = np.array([[-2.0, 1.0, 0.0],
                  [0.0, -2.0, 1.0],
                  [0.0, 0.0, -2.0]])

    # Step 3.1: Compute the cofactor matrix C_1.
    # For a 3x3 matrix, the adjugate (transpose of cofactor matrix) can be
    # calculated using the formula: adj(A) = det(A) * inv(A)
    # det(M) = (-2)^3 = -8
    det_M = np.linalg.det(M)
    # M is invertible since det(M) is not zero.
    inv_M = np.linalg.inv(M)
    Cofactor_T = det_M * inv_M
    C = Cofactor_T.T

    # Step 3.2: Compute the antisymmetric part of C.
    A_prime = (C - C.T) / 2

    # Step 3.3: Extract the tridiagonal part T.
    # This is our interpretation of "tridiagonal matrix of the Parlett-Reid decomposition"
    # for a skew-symmetric matrix with a zero diagonal.
    T = np.diag(np.diag(A_prime, -1), -1) + \
        np.diag(np.diag(A_prime, 0), 0) + \
        np.diag(np.diag(A_prime, 1), 1)

    # Step 3.4: Compute the square of T.
    T_squared = T @ T

    # Step 4: Calculate the largest Ky Fan norm (nuclear norm) of T_squared.
    # The singular values of a symmetric matrix are the absolute values of its eigenvalues.
    eigenvalues = np.linalg.eigvals(T_squared)
    singular_values = np.abs(eigenvalues)
    
    # Sort for deterministic output
    singular_values = np.sort(singular_values)[::-1]
    
    ky_fan_norm = np.sum(singular_values)
    
    # Print the equation with each number as requested.
    sv_strings = []
    for s in singular_values:
        # Format numbers to be clean (e.g., 2.0 -> 2)
        if np.isclose(s, round(s)):
            sv_strings.append(str(int(round(s))))
        else:
            sv_strings.append(f"{s:.4f}")

    final_sum_str = f"{ky_fan_norm:.4f}"
    if np.isclose(ky_fan_norm, round(ky_fan_norm)):
        final_sum_str = str(int(round(ky_fan_norm)))
        
    print(f"{' + '.join(sv_strings)} = {final_sum_str}")

solve_matrix_problem()