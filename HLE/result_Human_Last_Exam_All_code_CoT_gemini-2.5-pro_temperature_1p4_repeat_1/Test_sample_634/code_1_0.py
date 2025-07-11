import sympy

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two
    polynomial matrices P1 and P2.
    """
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])

    P2 = sympy.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    print("P1 =")
    sympy.pprint(P1)
    print("\nP2 =")
    sympy.pprint(P2)

    # Step 1: Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("\nInitial stacked matrix M:")
    sympy.pprint(M)

    # Step 2: Perform elementary row operations to reduce M.
    # We use rows of P2 (M[2,:] and M[3,:]) to eliminate elements in rows of P1 (M[0,:] and M[1,:]).
    # Notation: M.row_op(i, f) transforms row i to f(M.row(i), i)
    # Let's use R3 (M[3,:]) which has a 1 in the second column to clear that column.
    
    # R0 -> R0 + s*R3
    M.row_op(0, lambda val, j: sympy.expand(val + s * M[3, j]))
    
    # R1 -> R1 - s^2*R3
    M.row_op(1, lambda val, j: sympy.expand(val - s**2 * M[3, j]))

    print("\nAfter clearing elements using R3:")
    sympy.pprint(M)
    # M is now [[0, 0], [s**3 - 1, 0], [s, 0], [-s - 1, 1]]

    # Now use polynomial division on the first column for rows 1 and 2.
    # R1 -> R1 - s^2*R2
    M.row_op(1, lambda val, j: sympy.expand(val - s**2 * M[2, j]))
    
    print("\nAfter polynomial division (R1 -> R1 - s^2*R2):")
    sympy.pprint(M)
    # M is now [[0, 0], [-1, 0], [s, 0], [-s - 1, 1]]
    
    # Now use the new R1 = [-1, 0] to clear the first column of other rows.
    # R2 -> R2 + s*R1
    M.row_op(2, lambda val, j: sympy.expand(val + s * M[1, j]))
    
    # R3 -> R3 - (s+1)*R1
    M.row_op(3, lambda val, j: sympy.expand(val - (s + 1) * M[1, j]))
    
    print("\nAfter final clearing of the first column:")
    sympy.pprint(M)
    # M is now [[0, 0], [-1, 0], [0, 0], [0, 1]]

    # Step 3: Extract the non-zero rows to form the GCRD.
    # The non-zero rows are [-1, 0] and [0, 1].
    G_unnormalized = sympy.Matrix([M.row(1), M.row(3)])

    print("\nUnnormalized GCRD matrix G:")
    sympy.pprint(G_unnormalized)

    # Step 4: Normalize the GCRD.
    # The GCRD is unique up to left multiplication by a unimodular matrix.
    # We can multiply the first row by -1 to make the leading entry positive.
    # This corresponds to left-multiplying by [[-1, 0], [0, 1]].
    normalizer = sympy.Matrix([[-1, 0],[0, 1]])
    G_normalized = normalizer * G_unnormalized
    
    print("\nFinal normalized GCRD:")
    sympy.pprint(G_normalized)
    
    # Printing the final answer in the requested format
    print("\nThe final equation is:")
    print(f"GCRD = [[{G_normalized[0,0]}, {G_normalized[0,1]}], [{G_normalized[1,0]}, {G_normalized[1,1]}]]")


solve_gcrd()