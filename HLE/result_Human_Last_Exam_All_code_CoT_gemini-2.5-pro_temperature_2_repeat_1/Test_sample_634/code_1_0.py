import sympy

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices.
    """
    # Define the symbolic variable 's'
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

    # Step 1: Compute the transposes of P1 and P2
    P1_T = P1.transpose()
    P2_T = P2.transpose()

    # Step 2: Form the compound matrix M by stacking the transposes vertically
    M = P1_T.col_join(P2_T)

    # Step 3: Compute the row-reduced echelon form (RREF) of M.
    # The rref() method returns a tuple: (RREF_matrix, pivot_columns)
    # We only need the matrix part.
    rref_M = M.rref()[0]

    # The non-zero rows of the RREF form the GCLD of the transposed matrices.
    # Since P1 and P2 are 2x2, the GCRD will also be 2x2.
    # The GCLD will be the top 2 rows of the RREF matrix.
    # Note: Sympy's rref performs operations over the field of rational functions of 's',
    # which is the correct procedure for the polynomial Euclidean algorithm.
    num_cols = P1.shape[1]
    G_L = rref_M[:num_cols, :]

    # Step 4: The GCRD is the transpose of the GCLD found in the previous step.
    GCRD = G_L.transpose()

    # Print the resulting GCRD matrix element by element.
    print("The greatest common right divisor is a 2x2 matrix G(s) = [[g11, g12], [g21, g22]] where:")
    print(f"g11 = {GCRD[0, 0]}")
    print(f"g12 = {GCRD[0, 1]}")
    print(f"g21 = {GCRD[1, 0]}")
    print(f"g22 = {GCRD[1, 1]}")


solve_gcrd()