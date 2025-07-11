import sympy as sp

def solve_gcrd():
    """
    This function calculates the Greatest Common Right Divisor (GCRD)
    of two given polynomial matrices P1 and P2 using the extended
    Euclidean algorithm (row reduction method).
    """
    # Define the symbolic variable 's'
    s = sp.symbols('s')

    # Define the polynomial matrices P1 and P2
    P1 = sp.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])

    P2 = sp.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    print("Given polynomial matrices:")
    print("P1(s) =")
    sp.pprint(P1)
    print("\nP2(s) =")
    sp.pprint(P2)

    # --- Method Explanation ---
    print("\nTo find the Greatest Common Right Divisor (GCRD), we use the Euclidean algorithm for polynomial matrices.")
    print("We form a stacked matrix M(s) = [P1(s); P2(s)] and apply elementary polynomial row operations to reduce it to the form [G(s); 0].")
    print("The resulting matrix G(s) will be a GCRD of P1(s) and P2(s).\n")

    # Form the stacked matrix
    M = P1.col_join(P2)
    print("Initial stacked matrix M(s):")
    sp.pprint(M)
    
    # --- Row Operations ---
    print("\n--- Applying Elementary Row Operations ---")

    print("\nStep 1: Create a pivot with a constant leading term by adding row 3 to row 4 (R4 -> R4 + R3).")
    # M[3,:] corresponds to the 4th row (0-indexed)
    M[3,:] = M[3,:] + M[2,:]
    print("M(s) after Step 1:")
    sp.pprint(M)

    print("\nStep 2: Use the new row 4 (which is [-1, 1]) as a pivot to eliminate the first element in all other rows.")
    # R1 -> R1 + (s^2+s)*R4
    M[0,:] = M[0,:] + (s**2+s)*M[3,:]
    # R2 -> R2 - (s^2+1)*R4
    M[1,:] = M[1,:] - (s**2+1)*M[3,:]
    # R3 -> R3 + s*R4
    M[2,:] = M[2,:] + s*M[3,:]
    print("M(s) after Step 2:")
    sp.pprint(M)

    print("\nStep 3: Reorder the rows to get an upper triangular form for clarity.")
    # The current matrix M has been modified in place. Its rows are now:
    # [0, s**2], [0, -1], [0, s], [-1, 1]
    # We re-order them to place the non-zero pivot at the top: [-1, 1], [0, -1], [0, s], [0, s**2]
    M_reordered = sp.Matrix([ M.row(3), M.row(1), M.row(2), M.row(0) ])
    print("Reordered matrix:")
    sp.pprint(M_reordered)

    print("\nStep 4: Use the new upper block to eliminate the lower block, making the last two rows zero.")
    # R3 -> R3 + s*R2
    M_reordered[2,:] = M_reordered[2,:] + s*M_reordered[1,:]
    # R4 -> R4 + s^2*R2
    M_reordered[3,:] = M_reordered[3,:] + (s**2)*M_reordered[1,:]
    print("Matrix after final reduction:")
    sp.pprint(M_reordered)
    
    # --- Conclusion ---
    print("\nThe row operations have successfully reduced the stacked matrix to the form [G; 0].")
    # The GCRD is the top 2x2 non-zero block of the final matrix
    G = M_reordered[:2, :]
    print("The resulting GCRD matrix, G(s), is:")
    sp.pprint(G)

    print("\nThe elements of the final GCRD matrix G are:")
    print("G = [[{}, {}], [{}, {}]]".format(G[0,0], G[0,1], G[1,0], G[1,1]))

    det_G = G.det()
    print("\nThe determinant of G is {}, which is a non-zero constant.".format(det_G))
    print("A matrix with a constant non-zero determinant is called unimodular.")
    print("This means that P1 and P2 are right coprime, so their GCRD is any unimodular matrix.")
    print("A common choice for the GCRD in this case is the identity matrix I = [[1, 0], [0, 1]]. We will return the one computed from our procedure.")

# Run the solver
solve_gcrd()