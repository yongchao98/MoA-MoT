import numpy as np

def check_continuity_point(M):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.

    This is true if and only if the matrix is non-derogatory, which is equivalent
    to the degree of its minimal polynomial being n.
    
    The degree of the minimal polynomial is the dimension of the subspace spanned
    by the powers of M, i.e., dim(span{I, M, M^2, ...}).
    We can find this dimension by computing the rank of the matrix whose columns
    are the vectorized powers {I, M, M^2, ..., M^(n-1)}.
    """
    if M.shape[0] != M.shape[1]:
        print("Input is not a square matrix.")
        return

    n = M.shape[0]
    
    # Generate the powers of the matrix M^0, M^1, ..., M^(n-1)
    powers = []
    current_power = np.identity(n, dtype=M.dtype)
    for _ in range(n):
        powers.append(current_power.flatten())
        current_power = current_power @ M
        
    # Create a matrix from the vectorized powers
    # The vectors are columns, so we transpose the list of rows
    krylov_matrix = np.array(powers).T
    
    # The rank of this matrix is the degree of the minimal polynomial
    # Note: np.linalg.matrix_rank is robust for floating point numbers
    min_poly_degree = np.linalg.matrix_rank(krylov_matrix)
    
    # The matrix is a continuity point if the degree of the minimal poly is n
    is_continuity_point = (min_poly_degree == n)
    
    print(f"Testing matrix:\n{M}")
    print(f"Matrix dimension n: {n}")
    print(f"Degree of minimal polynomial: {min_poly_degree}")

    final_equation = f"{min_poly_degree} = {n}"
    
    if is_continuity_point:
        print(f"Is it a point of continuity? Yes, because {final_equation} is true.")
    else:
        print(f"Is it a point of continuity? No, because {final_equation} is false.")
    print("-" * 30)

if __name__ == '__main__':
    # Example 1: A non-derogatory matrix (a Jordan block).
    # This should be a point of continuity.
    # Its minimal polynomial is (x-2)^3, degree 3.
    M1 = np.array([[2, 1, 0],
                   [0, 2, 1],
                   [0, 0, 2]], dtype=float)
    check_continuity_point(M1)

    # Example 2: A derogatory matrix.
    # This should be a point of discontinuity.
    # Its minimal polynomial is (x-2)^2, degree 2 < 3.
    M2 = np.array([[2, 1, 0],
                   [0, 2, 0],
                   [0, 0, 2]], dtype=float)
    check_continuity_point(M2)
    
    # Example 3: A scalar matrix (which is derogatory).
    # This should be a point of discontinuity.
    # Its minimal polynomial is x-5, degree 1 < 4.
    M3 = np.identity(4) * 5
    check_continuity_point(M3)
