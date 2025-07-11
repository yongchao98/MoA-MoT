import sympy

def find_gcrd():
    """
    This function finds the greatest common right divisor (GCRD) of two
    polynomial matrices, P1 and P2, using the Euclidean algorithm based on
    elementary row operations.
    """
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s], [-s**2 - 1, s**2]])
    P2 = sympy.Matrix([[s, 0], [-s - 1, 1]])

    print("The given polynomial matrices are:")
    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)
    print("-" * 30)

    # Step 1: Form the augmented matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("Applying the Euclidean algorithm by performing row reduction on the composite matrix M = [P1; P2].")

    # Step 2: Apply a sequence of elementary row operations to reduce M.
    
    # Operations to zero out the first column below a certain pivot
    # R1 -> R1 - (s+1)*R3 and R2 -> R2 + s*R3
    M[0, :] = M[0, :] - (s + 1) * M[2, :]
    M[1, :] = M[1, :] + s * M[2, :]

    # Swap R1, R2 and make the leading element 1
    M.row_swap(0, 1)
    M[0, :] = -M[0, :]

    # Use the new R1 to eliminate other elements in the first column
    # R3 -> R3 - s*R1 and R4 -> R4 + (s+1)*R1
    M[2, :] = M[2, :] - s * M[0, :]
    M[3, :] = M[3, :] + (s + 1) * M[0, :]
    M[3, 1] = sympy.simplify(M[3, 1])

    # Operations to zero out the second column below a certain pivot
    # This part is analogous to the polynomial Euclidean algorithm
    # R3 -> R3 + s^2*R2 and R4 -> R4 - (s^2+s)*R2
    q, _ = sympy.div(sympy.poly(M[3, 1], s), sympy.poly(M[1, 1], s), domain='EX')
    M[2, :] = M[2, :] + (s**2) * M[1, :]
    M[3, :] = M[3, :] - q * M[1, :]

    # Final cleanup to get to row-echelon form
    M.row_swap(1, 3) # bring the row [0, 1] to the pivot position
    M[0, :] = M[0, :] + (s**2) * M[1, :] # R1 -> R1 + s^2*R2
    M[3, :] = M[3, :] + s * M[1, :] # R4 -> R4 + s*R2

    # Step 3: Extract the GCRD from the non-zero rows of the reduced matrix
    GCRD = M[:2, :]

    print("\nAfter reduction, the greatest common right divisor D(s) is:")
    print("D(s) = [[{}, {}],".format(GCRD[0,0], GCRD[0,1]))
    print("        [{}, {}]]".format(GCRD[1,0], GCRD[1,1]))

find_gcrd()