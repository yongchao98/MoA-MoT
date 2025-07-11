import sympy as sp

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices.
    """
    # Define the symbolic variable 's'
    s = sp.Symbol('s')

    # Define the polynomial matrices P1 and P2
    # P1 = [s^2 + s, -s;
    #       -s^2 - 1, s^2]
    # P2 = [s, 0;
    #       -s - 1, 1]

    # Form the stacked 4x2 matrix M = [P1; P2]
    M = sp.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2],
        [s, 0],
        [-s - 1, 1]
    ])

    # --- Perform Row Reduction (Euclidean Algorithm) ---

    # Swap row 0 and row 2 to get a lower-degree polynomial at the top for easier pivoting.
    M.row_swap(0, 2)
    # M is now [[s, 0], [-s**2 - 1, s**2], [s**2 + s, -s], [-s - 1, 1]]

    # Use the new row 0 to create zeros in the first column of the rows below it.
    M[1,:] = sp.simplify(M.row(1) + s * M.row(0))
    M[2,:] = sp.simplify(M.row(2) - (s + 1) * M.row(0))
    M[3,:] = sp.simplify(M.row(3) + 1 * M.row(0))
    # M is now [[s, 0], [-1, s**2], [0, -s], [-1, 1]]

    # Swap row 0 and row 1 to get -1 as the leading element, which is an excellent pivot.
    M.row_swap(0, 1)
    # M is now [[-1, s**2], [s, 0], [0, -s], [-1, 1]]

    # Use the new row 0 to eliminate the first element in other rows.
    M[1,:] = sp.simplify(M.row(1) + s * M.row(0))
    M[3,:] = sp.simplify(M.row(3) - 1 * M.row(0))
    # M is now [[-1, s**2], [0, s**3], [0, -s], [0, 1 - s**2]]

    # The first column is reduced. Now reduce the second column using the pivot row [0, -s].
    # Swap row 1 and row 2 to use the lower-degree pivot.
    M.row_swap(1, 2)
    # M is now [[-1, s**2], [0, -s], [0, s**3], [0, 1 - s**2]]
    
    # Use the new row 1 to eliminate elements in the second column of other rows.
    M[0,:] = sp.simplify(M.row(0) + s * M.row(1))
    M[2,:] = sp.simplify(M.row(2) + s**2 * M.row(1))
    M[3,:] = sp.simplify(M.row(3) - s * M.row(1))
    # M is now [[-1, 0], [0, -s], [0, 0], [0, 1]]

    # Use the last row [0, 1] to eliminate the remaining element in the second column.
    M[1,:] = sp.simplify(M.row(1) + s * M.row(3))
    # M is now [[-1, 0], [0, 0], [0, 0], [0, 1]]
    
    # --- Normalize the Result ---

    # The non-zero rows form the GCRD. Let's normalize them.
    # Multiply the first row by -1 to make the leading element 1.
    M[0,:] = -1 * M.row(0)
    # M is now [[1, 0], [0, 0], [0, 0], [0, 1]]

    # Form the GCRD matrix from the non-zero rows.
    gcrd_matrix = sp.Matrix([M.row(0), M.row(3)])

    # --- Print the Final Answer ---
    print("The greatest common right divisor (GCRD) is the matrix G:")
    
    # Print each element of the final matrix as requested.
    # The final GCRD is the identity matrix.
    print(f"G = [[{gcrd_matrix[0,0]}, {gcrd_matrix[0,1]}],")
    print(f"     [{gcrd_matrix[1,0]}, {gcrd_matrix[1,1]}]]")

solve_gcrd()