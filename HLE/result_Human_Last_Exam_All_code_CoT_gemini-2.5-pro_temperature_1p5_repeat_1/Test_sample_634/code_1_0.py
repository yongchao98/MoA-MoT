import sympy as sp

def solve_gcrd():
    """
    This function calculates the greatest common right divisor (GCRD) of two
    given polynomial matrices P1 and P2.
    """
    # Define 's' as a symbolic variable
    s = sp.Symbol('s')

    # Define the polynomial matrices P1 and P2 as given in the problem
    P1 = sp.Matrix([[s**2 + s, -s],
                    [-s**2 - 1, s**2]])

    P2 = sp.Matrix([[s, 0],
                    [-s - 1, 1]])

    # To find the GCRD, we stack P1 and P2 vertically to form a single matrix M
    M = sp.Matrix.vstack(P1, P2)

    # We then compute the row-reduced echelon form (RREF) of M.
    # The `rref()` method returns the RREF matrix and the indices of the pivot columns.
    M_rref, pivot_columns = M.rref()

    # The GCRD is formed by the non-zero rows of the RREF matrix.
    # The number of non-zero rows is equal to the number of pivot columns.
    num_gcrd_rows = len(pivot_columns)
    G = M_rref[:num_gcrd_rows, :]

    # Print the resulting GCRD matrix.
    # The output is formatted to clearly show the final matrix equation.
    print("The greatest common right divisor G(s) is:")
    
    # Build a string representation of the matrix for clean printing
    matrix_str = "G(s) = ["
    for i in range(G.rows):
        row_str = ", ".join(map(str, G.row(i)))
        matrix_str += f"[{row_str}]"
        if i < G.rows - 1:
            matrix_str += ",\n        "
    matrix_str += "]"
    
    print(matrix_str)

solve_gcrd()