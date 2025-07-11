import sympy

def solve():
    """
    This function calculates the greatest common right divisor (GCRD) of two
    polynomial matrices P1 and P2 using the row reduction method.
    """
    # Define the symbol for the polynomial variable
    s = sympy.Symbol('s')

    # Define the two polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])
    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # To find the GCRD, we stack the matrices and perform elementary row operations
    # to find a basis for the row space.
    # Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)

    # Create a copy of the matrix to perform row reductions
    M_reduced = M.copy()

    # Perform a series of row operations to bring the matrix to a row-echelon form.
    # This is analogous to Gaussian elimination for polynomial matrices.
    # The exact sequence of operations can vary, but the result for the GCRD basis is the same.
    # Step 1: Use row 2 (index 2 of M) to simplify other rows.
    M_reduced.row_op(0, lambda val, j: sympy.expand(val - (s + 1) * M_reduced[2, j]))
    M_reduced.row_op(1, lambda val, j: sympy.expand(val + s * M_reduced[2, j]))
    M_reduced.row_op(3, lambda val, j: sympy.expand(val + M_reduced[2, j]))
    
    # Step 2: Use the new row 3 (which is now [-1, 1]) as a pivot.
    M_reduced.row_swap(0, 3) # Swap row 0 and 3
    M_reduced.row_op(0, lambda val, j: -val) # Multiply the new row 0 by -1 to get a leading 1
    
    # Step 3: Eliminate other elements in the first column using the new pivot row 0.
    M_reduced.row_op(1, lambda val, j: sympy.expand(val + M_reduced[0, j]))
    M_reduced.row_op(2, lambda val, j: sympy.expand(val - s * M_reduced[0, j]))
    
    # Step 4: Reduce the remaining sub-matrix in the second column.
    M_reduced.row_op(3, lambda val, j: sympy.expand(val + M_reduced[2, j]))
    M_reduced.row_op(1, lambda val, j: sympy.expand(val - s * M_reduced[2, j]))
    M_reduced.row_op(2, lambda val, j: sympy.expand(val + s * M_reduced[1, j]))

    # After reduction, the matrix has two non-zero rows which form a GCRD.
    # These rows are [1, -1] and [0, -1]. We extract them.
    GCRD = M_reduced.minor_submatrix(0, 0, 2, 2)

    # To get the simplest, canonical form (row-echelon form), we do final reductions.
    GCRD.row_op(1, lambda val, j: -val)  # R1 -> -R1, GCRD becomes [[1, -1], [0, 1]]
    GCRD.row_op(0, lambda val, j: sympy.expand(val + GCRD[1, j])) # R0 -> R0 + R1

    # The resulting GCRD is the identity matrix.
    final_gcrd = GCRD

    print("The greatest common right divisor (GCRD) of P1 and P2 is:")

    # Print each element of the final matrix as requested.
    row1 = [final_gcrd[0, 0], final_gcrd[0, 1]]
    row2 = [final_gcrd[1, 0], final_gcrd[1, 1]]
    print(f"[{row1[0]}, {row1[1]}]")
    print(f"[{row2[0]}, {row2[1]}]")

solve()