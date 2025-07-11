import sympy

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices
    P1 and P2 by performing elementary row operations on the stacked matrix [P1; P2].
    """
    # Define the symbolic variable 's'
    s = sympy.symbols('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])
    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    print("Given matrices:")
    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)

    # Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("\nStep 1: Form the stacked matrix M = [P1; P2]")
    sympy.pprint(M)

    # To simplify, we want to use lower-degree rows to eliminate higher-degree rows.
    # Swap rows to bring the lower-degree rows from P2 to the top.
    print("\nStep 2: Swap rows to bring lower-degree polynomials to the top for pivoting.")
    M.row_swap(0, 2)
    M.row_swap(1, 3)
    print("M after swapping rows:")
    sympy.pprint(M)

    # Now M = [[s, 0], [-s - 1, 1], [s^2 + s, -s], [-s^2 - 1, s^2]]
    # Use the first row [s, 0] to eliminate entries in the first column.
    print("\nStep 3: Eliminate elements in the first column.")
    # R3 := R3 - (s+1)*R1
    M[2, :] = M.row(2) - (s + 1) * M.row(0)
    # R4 := R4 + s*R1
    M[3, :] = M.row(3) + s * M.row(0)
    print("M after eliminating s^2 terms:")
    sympy.pprint(M)

    # Now M = [[s, 0], [-s - 1, 1], [0, -s], [-1, s^2]]
    # A '1' or '-1' is an excellent pivot. Let's use the new R4 [-1, s^2].
    print("\nStep 4: Use the row [-1, s^2] to clear the rest of the first column.")
    # R1 := R1 + s*R4
    M[0, :] = M.row(0) + s * M.row(3)
    # R2 := R2 - (s+1)*R4
    M[1, :] = M.row(1) - (s + 1) * M.row(3)
    print("M after clearing the first column:")
    sympy.pprint(M)

    # Now M = [[0, s^3], [0, 1 - s^2 - s^3], [0, -s], [-1, s^2]]
    # Let's clean up by moving the non-zero row in the first column to the top.
    print("\nStep 5: Rearrange rows and normalize the pivot row.")
    M.row_swap(0, 3)
    M[0, :] = -M.row(0) # R1 := -1 * R1
    print("M after rearranging:")
    sympy.pprint(M)

    # Now M = [[1, -s^2], [0, 1 - s^2 - s^3], [0, -s], [0, s^3]]
    # The first column is done. Now we work on the submatrix in the second column.
    # We use the Euclidean algorithm on the polynomials in the second column.
    # The GCD of (1 - s^2 - s^3), (-s), and (s^3) is 1. We can get a row [0, 1].
    print("\nStep 6: Reduce the second column using row operations.")
    # First, use R3 = [0, -s] to simplify R4 = [0, s^3]
    # R4 := R4 + s^2*R3
    M[3, :] = M.row(3) + (s**2) * M.row(2)
    # Now use R3 to simplify R2 = [0, 1 - s^2 - s^3]
    # From polynomial division, (1 - s^2 - s^3) = (s^2 + s)*(-s) + 1.
    # So, 1 = (1 - s^2 - s^3) - (s^2 + s)*(-s).
    # R2_new = R2 - (s^2+s)*R3
    M[1, :] = M.row(1) - (s**2 + s) * M.row(2)
    print("M after reducing with [0, -s]:")
    sympy.pprint(M)

    # Now M = [[1, -s^2], [0, 1], [0, -s], [0, 0]]
    # We have a row [0, 1], which is an excellent pivot for the second column.
    print("\nStep 7: Use the row [0, 1] to clear the rest of the second column.")
    # R1 := R1 + s^2 * R2
    M[0, :] = M.row(0) + (s**2) * M.row(1)
    # R3 := R3 + s * R2
    M[2, :] = M.row(2) + s * M.row(1)
    print("M after final reduction:")
    sympy.pprint(M)

    # The resulting matrix is [[1, 0], [0, 1], [0, 0], [0, 0]].
    # The non-zero part is the GCRD.
    GCRD = M[:2, :]
    print("\nThe row reduction is complete. The non-zero part is the GCRD.")
    print("The Greatest Common Right Divisor is:")
    
    # Print the final result element by element as requested
    final_eq = f"[[{GCRD[0,0]}, {GCRD[0,1]}],\n [{GCRD[1,0]}, {GCRD[1,1]}]]"
    print(final_eq)

solve_gcrd()