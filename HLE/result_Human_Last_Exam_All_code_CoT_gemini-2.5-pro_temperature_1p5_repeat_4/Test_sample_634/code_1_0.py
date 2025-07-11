import sympy as sp

def find_gcrd(P1, P2):
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices.
    """
    s = sp.Symbol('s')

    # Step 1 & 2: Transpose the input matrices
    P1_T = P1.transpose()
    P2_T = P2.transpose()
    print("P1 Transposed:")
    sp.pprint(P1_T)
    print("\nP2 Transposed:")
    sp.pprint(P2_T)
    print("-" * 30)
    
    # Step 3: Find the GCLD of the transposed matrices using row reduction
    # Form the stacked matrix M = [P1^T; P2^T]
    M = P1_T.col_join(P2_T)
    print("Stacked matrix M = [P1^T; P2^T]:")
    sp.pprint(M)
    print("-" * 30)

    # Make a copy to perform row operations on
    M_reduced = M.copy()

    # The last row R4 is [0, 1]. Use it to clear the second column of other rows.
    # Note: Sympy rows are 0-indexed, so R1 is row 0, R2 is row 1, etc.
    # R1 -> R1 + (s^2 + 1) * R4
    M_reduced.row_op(0, lambda val, j: sp.expand(val + (s**2 + 1) * M_reduced[3, j]))
    # R2 -> R2 - s^2 * R4
    M_reduced.row_op(1, lambda val, j: sp.expand(val - s**2 * M_reduced[3, j]))
    # R3 -> R3 + (s + 1) * R4
    M_reduced.row_op(2, lambda val, j: sp.expand(val + (s + 1) * M_reduced[3, j]))

    print("After clearing the second column using R4:")
    sp.pprint(M_reduced)
    print("-" * 30)

    # Now the matrix is:
    # [s^2 + s, 0]
    # [   -s,   0]
    # [    s,   0]
    # [    0,   1]
    
    # Reduce the first three rows. We can eliminate R1 and R2 using R3.
    # The GCD of s^2+s, -s, s is s.
    # R1 -> R1 - (s+1) * R3
    M_reduced.row_op(0, lambda val, j: sp.expand(val - (s + 1) * M_reduced[2, j]))
    # R2 -> R2 + R3
    M_reduced.row_op(1, lambda val, j: sp.expand(val + M_reduced[2, j]))
    
    print("After reducing rows with common factor:")
    sp.pprint(M_reduced)
    print("-" * 30)

    # The final non-zero rows form the GCLD of the transposed matrices.
    # These rows are [s, 0] and [0, 1].
    # We can form the GCLD matrix G_L in row-echelon form.
    GCLD_T = sp.Matrix([M_reduced.row(2), M_reduced.row(3)])
    
    print("GCLD of transposed matrices is:")
    sp.pprint(GCLD_T)
    print("-" * 30)
    
    # Step 4: The GCRD is the transpose of GCLD_T
    GCRD = GCLD_T.transpose()
    return GCRD

if __name__ == '__main__':
    s = sp.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sp.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])
    P2 = sp.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    # Calculate the GCRD
    final_gcrd = find_gcrd(P1, P2)

    print("The final greatest common right divisor (GCRD) matrix is:")
    sp.pprint(final_gcrd)

    # Output each element of the final GCRD matrix
    print("\nThe elements of the GCRD matrix are:")
    for r in range(final_gcrd.rows):
        for c in range(final_gcrd.cols):
            print(f"G[{r},{c}] = {final_gcrd[r, c]}")
