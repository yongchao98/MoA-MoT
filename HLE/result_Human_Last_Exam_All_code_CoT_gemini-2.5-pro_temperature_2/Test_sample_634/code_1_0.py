import sympy as sp

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices
    P1 and P2 using the elementary row operations method.
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

    print("Original Matrix P1:")
    sp.pretty_print(P1)
    print("\nOriginal Matrix P2:")
    sp.pretty_print(P2)

    # Step 1: Form the stacked matrix M
    M = P1.col_join(P2)
    print("\nInitial stacked matrix M:")
    sp.pretty_print(M)
    
    # Step 2: Perform elementary row operations
    
    # Let the rows of M be R0, R1, R2, R3
    # M = [ R0 ]
    #     [ R1 ]
    #     [ R2 ]
    #     [ R3 ]
    
    # Reduce rows with the highest degree (deg=2) using rows with lower degree (deg=1)
    # R0_new = R0 - (s+1)*R2
    M[0, :] = M.row(0) - (s + 1) * M.row(2)
    # R1_new = R1 + s*R2
    M[1, :] = M.row(1) + s * M.row(2)
    print("\nAfter first reduction step (reducing degrees of top rows):")
    sp.pretty_print(M)
    # M is now:
    # [[   0,   -s],
    #  [  -1,  s**2],
    #  [   s,     0],
    #  [-s-1,     1]]

    # Reduce the new R1 (now at index 1), which has the highest degree (2).
    # R1_new = R1 + s * R0
    M[1, :] = M.row(1) + s * M.row(0)
    print("\nAfter second reduction step (reducing row 1):")
    sp.pretty_print(M)
    # M is now:
    # [[   0, -s],
    #  [  -1,  0],
    #  [   s,  0],
    #  [-s-1,  1]]

    # Reduce rows 2 and 3 using the new row 1 (at index 1), which is now degree 0.
    # R2_new = R2 + s * R1
    M[2, :] = M.row(2) + s * M.row(1)
    # R3_new = R3 - (s+1) * R1
    M[3, :] = M.row(3) - (s + 1) * M.row(1)
    print("\nAfter third reduction step (reducing rows 2 and 3):")
    sp.pretty_print(M)
    # M is now:
    # [[ 0, -s],
    #  [-1,  0],
    #  [ 0,  0],
    #  [ 0,  1]]

    # Reduce row 0 using the new row 3.
    # R0_new = R0 + s * R3
    M[0, :] = M.row(0) + s * M.row(3)
    print("\nAfter final reduction step (reducing row 0):")
    sp.pretty_print(M)
    # M is now:
    # [[0,  0],
    #  [-1, 0],
    #  [0,  0],
    #  [0,  1]]

    # Step 3: Extract the GCRD
    # The non-zero rows form the GCRD.
    # The non-zero rows are [-1, 0] and [0, 1].
    G_unnormalized = sp.Matrix([M.row(1), M.row(3)])
    
    # Step 4: Normalize the result
    # We can multiply rows by -1 to make the diagonal elements positive (monic).
    G_unnormalized[0,:] = -1 * G_unnormalized.row(0)
    G = sp.Matrix([G_unnormalized.row(0), G_unnormalized.row(1)])

    # Rearranging for canonical form I
    if G.row(0) == sp.Matrix([[1,0]]) and G.row(1) == sp.Matrix([[0,1]]):
        GCRD = sp.Matrix([[1, 0], [0, 1]])
    else: # This case is for ordering rows if needed, but here it's already correct.
        GCRD = G

    print("\nThe GCRD matrix G(s) is:")
    sp.pretty_print(GCRD)
    
    print("\nThe determinant of the GCRD is:", GCRD.det())
    print("Since the determinant is a non-zero constant, the matrices P1 and P2 are right coprime.")
    
    # As requested, output each number in the final equation for the GCRD matrix
    print("\nFinal equation for the GCRD matrix:")
    for i in range(GCRD.rows):
        for j in range(GCRD.cols):
            print(f"G[{i},{j}] = {GCRD[i,j]}")
            
    return GCRD

G = solve_gcrd()
<<<[[1, 0], [0, 1]]>>>