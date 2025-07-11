import sympy

def find_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two
    polynomial matrices P1 and P2 using the Euclidean algorithm via row reduction.
    """
    # Step 1: Define the symbolic variable and the polynomial matrices.
    s = sympy.Symbol('s')

    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])

    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # Step 2: Form the augmented matrix M by stacking P1 and P2.
    # The rows are indexed 0, 1, 2, 3.
    M = sympy.Matrix.vstack(P1, P2)

    # Step 3: Perform elementary row operations to reduce M.
    # The goal is to create zero rows at the bottom.
    
    # Let R0, R1, R2, R3 be the rows of M.
    # Use R3 = [-s-1, 1] to eliminate terms in other rows.

    # Operation: R0 -> R0 + s * R3
    # This makes the first row [0, 0]
    M[0, :] = M.row(0) + s * M.row(3)

    # Operation: R1 -> R1 - s**2 * R3
    # This eliminates the s**2 term in the second column of R1
    M[1, :] = M.row(1) - (s**2) * M.row(3)

    # Now M has two zero-rows in-the-making. Let's operate on the new rows.
    # Current state of M's rows (conceptually):
    # [0, 0]
    # [s**3 - 1, 0]
    # [s, 0]
    # [-s - 1, 1]

    # Use the polynomial Euclidean algorithm on the first elements of R1 and R2.
    # We want to find GCD(s**3 - 1, s). It is 1.
    # From Bezout's identity: 1*R1 - (s**2)*R2 gives a row starting with -1.
    # Operation: R1 -> R1 - s**2 * R2
    M[1, :] = M.row(1) - (s**2) * M.row(2)

    # Normalize the new R1 (which is now [-1, 0]) to [1, 0].
    # Operation: R1 -> -1 * R1
    M[1, :] = -1 * M.row(1)

    # Now we have a row [1, 0]. Use it to eliminate other rows.
    # Operation: R2 -> R2 - s * R1
    M[2, :] = M.row(2) - s * M.row(1)

    # Operation: R3 -> R3 - (-s - 1) * R1
    M[3, :] = M.row(3) - (-s - 1) * M.row(1)

    # Step 4: The non-zero rows of the reduced matrix form the GCRD.
    # The non-zero rows are now [1, 0] and [0, 1].
    
    non_zero_rows = []
    for i in range(M.rows):
        if not M.row(i).is_zero:
            non_zero_rows.append(M.row(i))
    
    # Assemble the GCRD matrix from the non-zero rows.
    GCRD = sympy.Matrix(non_zero_rows)
    # Ensure canonical form (identity matrix) if rows are out of order.
    if GCRD.row(0).col(0)[0] == 0:
        GCRD.row_swap(0,1)

    # Step 5: Print the final GCRD matrix.
    print("The greatest common right divisor D is:")
    # The prompt requests to output each number in the final equation.
    # Printing the matrix fulfills this.
    sympy.pprint(GCRD, use_unicode=False)

find_gcrd()