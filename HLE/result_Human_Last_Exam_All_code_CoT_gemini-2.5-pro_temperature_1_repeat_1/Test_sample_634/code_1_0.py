import sympy

def solve_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices.
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])
    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)
    print("-" * 30)

    # 1. Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("Step 1: Form the stacked matrix M = [P1; P2]")
    sympy.pprint(M)
    print("-" * 30)

    # 2. Perform elementary row operations to reduce M
    # Let the rows be R1, R2, R3, R4
    # M =
    # [R1: s**2 + s,    -s  ]
    # [R2: -s**2 - 1,   s**2]
    # [R3:     s,        0  ]
    # [R4:   -s - 1,     1  ]
    
    print("Step 2: Reduce M using elementary row operations.")
    
    # R1 -> R1 - (s+1)*R3
    # R2 -> R2 + s*R3
    # R4 -> R4 + R3
    M[0, :] = M[0, :] - (s + 1) * M[2, :]
    M[1, :] = M[1, :] + s * M[2, :]
    M[3, :] = M[3, :] + M[2, :]
    # print("After first set of operations:")
    # sympy.pprint(M)

    # Swap R2 and R4 to get a simpler pivot [-1, 1]
    rows = M.tolist()
    rows[1], rows[3] = rows[3], rows[1]
    M = sympy.Matrix(rows)
    # print("After swapping rows:")
    # sympy.pprint(M)

    # R3 -> R3 + s*R2 (new R2)
    # R4 -> R4 - R2 (new R2)
    M[2, :] = M[2, :] + s * M[1, :]
    M[3, :] = M[3, :] - M[1, :]
    # print("After second set of operations:")
    # sympy.pprint(M)

    # R1 -> R1 + R3
    # R4 -> R4 - s*R3
    M[0, :] = M[0, :] + M[2, :]
    M[3, :] = M[3, :] - s * M[2, :]
    # print("After third set of operations:")
    # sympy.pprint(M)
    
    # R2 -> R2 + R4 (new R4)
    # R3 -> R3 + s*R4 (new R4)
    M[1, :] = M[1, :] + M[3, :]
    M[2, :] = M[2, :] + s * M[3, :]
    
    print("\nFinal reduced matrix:")
    sympy.pprint(M)
    print("-" * 30)
    
    # 3. The non-zero rows of the reduced matrix form a GCRD.
    # The non-zero rows are [-1, 0] and [0, -1].
    G_derived = sympy.Matrix([M[1,:], M[3,:]])
    
    print("Step 3: A GCRD is formed by the non-zero rows of the reduced matrix.")
    print("Derived GCRD =")
    sympy.pprint(G_derived)
    print("-" * 30)

    # 4. Normalize the GCRD.
    # The derived GCRD is [-1, 0; 0, -1]. Its determinant is 1, so it is unimodular.
    # Any unimodular matrix is a valid GCRD for right-coprime matrices.
    # The simplest and most common representation is the identity matrix.
    # We can obtain it by multiplying by the unimodular matrix [[-1, 0], [0, -1]].
    U = sympy.Matrix([[-1, 0], [0, -1]])
    GCRD = U * G_derived
    
    print("Step 4: Normalize the GCRD to its simplest form.")
    print("The final GCRD is the identity matrix, which means P1 and P2 are right coprime.")
    
    # Print the final result, showing each element as requested.
    print("\nThe greatest common right divisor is:")
    print(f"[[{GCRD[0,0]}, {GCRD[0,1]}],")
    print(f" [{GCRD[1,0]}, {GCRD[1,1]}]]")

solve_gcrd()