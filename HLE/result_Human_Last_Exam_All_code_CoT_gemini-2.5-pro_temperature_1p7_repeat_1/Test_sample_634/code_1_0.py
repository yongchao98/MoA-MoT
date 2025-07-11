import sympy

def find_gcrd():
    """
    This function finds the Greatest Common Right Divisor (GCRD) of two
    given polynomial matrices P1 and P2 using the Euclidean algorithm method.
    """
    # Set up the symbolic variable and matrices
    s = sympy.symbols('s')

    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])

    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # Explain the initial setup
    print("The initial matrices are P1 and P2:")
    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)

    # Form the composite matrix M
    M = P1.col_join(P2)
    print("\nWe form the composite matrix M = [P1; P2]:")
    sympy.pprint(M)

    # --- Start of Elementary Row Operations ---

    # The 4th row has a constant '1' at M[3, 1], which makes a good pivot.
    print("\nStep 1: Use R4 to eliminate elements in the second column of R1 and R2.")
    print("Performing R1 <- R1 + s*R4 and R2 <- R2 - s^2*R4")
    
    # R1 <- R1 + s*R4
    M[0, :] = sympy.expand(M[0, :] + s * M[3, :])
    # R2 <- R2 - s^2*R4
    M[1, :] = sympy.expand(M[1, :] - s**2 * M[3, :])

    print("The matrix becomes:")
    sympy.pprint(M)

    # After Step 1, the matrix has zero rows and rows with polynomials only in the first column.
    # We now apply the Euclidean algorithm to the polynomials in the first column of the non-zero rows.
    # The relevant polynomials are from R2 (s**3 - 1) and R3 (s).
    print("\nStep 2: Apply Euclidean algorithm on the first column of R2 and R3.")
    print("Performing R2 <- R2 - s^2*R3, since (s^3 - 1) = s^2*(s) - 1.")

    M[1, :] = sympy.expand(M[1, :] - s**2 * M[2, :])
    
    print("The matrix becomes:")
    sympy.pprint(M)

    # Now we have a constant (-1) in R2, which is an excellent pivot to eliminate other elements.
    print("\nStep 3: Use the new R2 to eliminate other elements in the first column.")
    print("Performing R3 <- R3 + s*R2 and R4 <- R4 - (s+1)*R2")

    # R3 <- R3 + s*R2
    M[2, :] = sympy.expand(M[2, :] + s * M[1, :])
    # R4 <- R4 - (s+1)*R2
    M[3, :] = sympy.expand(M[3, :] - (s + 1) * M[1, :])
    
    print("The matrix becomes:")
    sympy.pprint(M)

    # The reduction is complete. The non-zero rows form the GCRD.
    # These rows are R2 = [-1, 0] and R4 = [0, 1].
    print("\nThe non-zero rows form the GCRD. We arrange them into a matrix.")
    G = sympy.Matrix([M[1,:], M[3,:]])

    # For a canonical form, we make the leading diagonal elements positive.
    if G[0, 0] == -1:
        G[0, :] = -G[0, :]

    print("After normalizing to make the leading coefficient positive, the GCRD is:")
    sympy.pprint(G)

    print("\nFinal Answer: The greatest common right divisor is the matrix G(s) shown above.")
    print("Each number in the final matrix is:")
    print(f"Row 1, Col 1: {G[0,0]}")
    print(f"Row 1, Col 2: {G[0,1]}")
    print(f"Row 2, Col 1: {G[1,0]}")
    print(f"Row 2, Col 2: {G[1,1]}")

if __name__ == '__main__':
    find_gcrd()