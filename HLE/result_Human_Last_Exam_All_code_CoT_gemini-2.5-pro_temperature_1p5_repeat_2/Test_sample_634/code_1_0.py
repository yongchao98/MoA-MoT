import sympy as sp

def find_gcrd():
    """
    This function finds the greatest common right divisor (GCRD) of two
    polynomial matrices P1 and P2 using the extended Euclidean algorithm.
    """
    # Define the symbolic variable 's'
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

    print("Original Matrix P1:")
    sp.pprint(P1)
    print("\nOriginal Matrix P2:")
    sp.pprint(P2)

    # Step 1: Form the augmented matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("\nAugmented Matrix M = [P1; P2]:")
    sp.pprint(M)

    # Step 2: Perform elementary row operations to reduce M.
    # The goal is to get a matrix of the form [G; 0].
    # To simplify, we'll swap rows to get the simpler matrix P2 on top.
    # This is equivalent to finding GCRD(P2, P1), which is the same as GCRD(P1, P2).
    M_reduced = M.copy()
    M_reduced.row_swap(0, 2)
    M_reduced.row_swap(1, 3)
    # M_reduced is now [P2; P1]

    # R3_new = R3 - (s+1)*R1
    # R4_new = R4 + s*R1
    # Note: Matrix rows are 0-indexed.
    M_reduced[2, :] = M_reduced.row(2) - (s + 1) * M_reduced.row(0)
    M_reduced[3, :] = M_reduced.row(3) + s * M_reduced.row(0)
    # This zeroes out the top-left element of the bottom block.

    # R2_new = R2 + R1
    M_reduced[1, :] = M_reduced.row(1) + M_reduced.row(0)
    # This simplifies row 1 to [-1, 1]

    # Now we have a very simple row [-1, 1]. Let's move it to the top.
    M_reduced.row_swap(0, 1)
    M_reduced[0, :] = -1 * M_reduced.row(0) # Make leading element 1.
    # Top row is now [1, -1]

    # Use the new top row [1, -1] to eliminate the first column of other rows.
    M_reduced[1, :] = M_reduced.row(1) - s * M_reduced.row(0)
    M_reduced[3, :] = M_reduced.row(3) + M_reduced.row(0)
    # First column is now [1, 0, 0, 0]^T

    # The matrix is now:
    # [ 1,   -1   ]
    # [ 0,    s   ]
    # [ 0,   -s   ]
    # [ 0, s**2-1 ]
    # Now, reduce the remaining rows.
    
    # R3_new = R3 + R2
    M_reduced[2, :] = M_reduced.row(2) + M_reduced.row(1)
    # Row 3 is now [0, 0]

    # R4_new = R4 - s*R2 (Polynomial division: (s^2-1) = s*(s) - 1)
    M_reduced[3, :] = M_reduced.row(3) - s * M_reduced.row(1)
    # Row 4 is now [0, -1]
    
    # The matrix is now:
    # [ 1,  -1 ]
    # [ 0,   s ]
    # [ 0,   0 ]
    # [ 0,  -1 ]
    # We have a new simple row [0, -1]. Use it to simplify further.

    # R2_new = R2 + s*R4
    M_reduced[1, :] = M_reduced.row(1) + s * M_reduced.row(3)
    # Row 2 is now [0, 0]

    # R1_new = R1 - R4
    M_reduced[0, :] = M_reduced.row(0) - M_reduced.row(3)
    # Row 1 is now [1, 0]

    # The matrix is now:
    # [ 1,  0 ]
    # [ 0,  0 ]
    # [ 0,  0 ]
    # [ 0, -1 ]
    # The non-zero rows form the GCRD.
    
    # Let's form the GCRD from the non-zero rows [1, 0] and [0, -1].
    # To put it in a canonical form (identity matrix), we multiply the last row by -1.
    G = sp.Matrix([
        M_reduced.row(0),
        -1 * M_reduced.row(3)
    ])

    print("\nFinal GCRD Matrix:")
    # Print each element as requested
    print(f"[{G[0,0]}, {G[0,1]}]")
    print(f"[{G[1,0]}, {G[1,1]}]")

if __name__ == '__main__':
    find_gcrd()