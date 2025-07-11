import sympy

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices.
    """
    # Define the symbolic variable 's'
    s = sympy.symbols('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])

    P2 = sympy.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    print("Given matrices:")
    print("P1 =")
    sympy.pprint(P1)
    print("\nP2 =")
    sympy.pprint(P2)
    
    # Step 1: Form the augmented matrix M = [P1; P2]
    M = P1.vstack(P2)
    print("\nStep 1: Form the augmented matrix M by stacking P1 and P2.")
    print("M =")
    sympy.pprint(M)

    # Let's denote the rows as R0, R1, R2, R3
    # M = [ R0 ]
    #     [ R1 ]
    #     [ R2 ]
    #     [ R3 ]
    # R3 is [-s-1, 1]. The element '1' is a convenient pivot.

    # Step 2: Use R3 to eliminate the second column entries in other rows.
    # R0_new = R0 + s * R3
    # R1_new = R1 - s^2 * R3
    print("\nStep 2: Use Row 3 ([-s - 1, 1]) to create zeros in the second column.")
    M_step2 = M.copy()
    M_step2[0, :] = M_step2[0, :] + s * M_step2[3, :]
    M_step2[1, :] = M_step2[1, :] - (s**2) * M_step2[3, :]
    print("Matrix after using R3 as a pivot:")
    sympy.pprint(M_step2)

    # The resulting matrix is:
    # [ 0,         0       ]
    # [ s**3 - 1,  0       ]
    # [ s,         0       ]
    # [ -s - 1,    1       ]

    # Step 3: Now we focus on the first column. Use R2 ([s, 0]) to simplify R1 ([s**3-1, 0]).
    # R1_new = R1 - s^2 * R2
    print("\nStep 3: Use the new Row 2 ([s, 0]) to simplify the new Row 1 ([s**3 - 1, 0]).")
    M_step3 = M_step2.copy()
    M_step3[1, :] = M_step3[1, :] - (s**2) * M_step3[2, :]
    print("Matrix after simplifying the first column:")
    sympy.pprint(M_step3)
    
    # The resulting matrix is:
    # [ 0,      0 ]
    # [ -1,     0 ]
    # [ s,      0 ]
    # [ -s - 1, 1 ]

    # Step 4: Use the new simpler R1 ([-1, 0]) to eliminate other rows.
    # R2_new = R2 + s * R1
    # R3_new = R3 - (s+1) * R1
    print("\nStep 4: Use the very simple Row 1 ([-1, 0]) to simplify Rows 2 and 3.")
    M_step4 = M_step3.copy()
    M_step4[2, :] = M_step4[2, :] + s * M_step4[1, :]
    M_step4[3, :] = M_step4[3, :] - (s + 1) * M_step4[1, :]
    print("Matrix after final reductions:")
    sympy.pprint(M_step4)
    
    # The resulting matrix is:
    # [ 0,  0 ]
    # [ -1, 0 ]
    # [ 0,  0 ]
    # [ 0,  1 ]
    
    # Step 5: Collect the non-zero rows and normalize them.
    # The non-zero rows are [-1, 0] and [0, 1].
    # We can multiply [-1, 0] by -1 (a unimodular operation) to get [1, 0].
    # Arranging them gives the GCRD.
    GCRD = sympy.Matrix([
        [-1 * M_step4[1, :]],
        [M_step4[3, :]]
    ])
    
    print("\nStep 5: The non-zero rows form the GCRD. After normalization (multiplying the first row by -1), we get:")
    print("GCRD =")
    sympy.pprint(GCRD)
    
    print("\nThe greatest common right divisor is the 2x2 identity matrix.")
    print("The final matrix elements are:")
    print(f"Row 1: {GCRD[0,0]}, {GCRD[0,1]}")
    print(f"Row 2: {GCRD[1,0]}, {GCRD[1,1]}")

solve_gcrd()