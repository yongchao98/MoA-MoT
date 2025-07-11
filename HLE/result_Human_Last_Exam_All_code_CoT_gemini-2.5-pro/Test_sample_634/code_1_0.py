import sympy

def solve():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices.
    """
    # Set up the symbolic variable and matrices
    s = sympy.Symbol('s')
    P1 = sympy.Matrix([[s**2 + s, -s], [-s**2 - 1, s**2]])
    P2 = sympy.Matrix([[s, 0], [-s - 1, 1]])

    # The GCRD is found by forming an augmented matrix M = [P1; P2] and
    # reducing it to the form [G; 0] using elementary row operations.

    # Form the augmented matrix M
    M = P1.col_join(P2)

    # Step 1: Use row 3 (original M[2,:]) to reduce rows 1 and 2.
    # R1_new = R1 - (s+1)*R3
    M.row_op(0, lambda val, j: sympy.expand(val - (s + 1) * M[2, j]))
    # R2_new = R2 + s*R3
    M.row_op(1, lambda val, j: sympy.expand(val + s * M[2, j]))

    # Step 2: Use row 3 to reduce row 4.
    # R4_new = R4 + R3
    M.row_op(3, lambda val, j: sympy.expand(val + M[2, j]))
    
    # At this point, the matrix is:
    # [  0   ,   -s  ]
    # [ -1   ,  s**2 ]
    # [  s   ,    0  ]
    # [ -1   ,    1  ]

    # Step 3: Reorder rows to bring simpler rows (with constants) to the top for pivoting.
    M = sympy.Matrix([M.row(3), M.row(1), M.row(2), M.row(0)])
    
    # At this point, the matrix is:
    # [ -1   ,    1  ]
    # [ -1   ,  s**2 ]
    # [  s   ,    0  ]
    # [  0   ,   -s  ]

    # Step 4: Use the new row 1 ([-1, 1]) to eliminate elements in the first column.
    # R2_new = R2 - R1
    M.row_op(1, lambda val, j: sympy.expand(val - M[0, j]))
    # R3_new = R3 + s*R1
    M.row_op(2, lambda val, j: sympy.expand(val + s * M[0, j]))

    # At this point, the matrix is:
    # [ -1   ,     1    ]
    # [  0   , s**2 - 1 ]
    # [  0   ,     s    ]
    # [  0   ,    -s    ]

    # Step 5: Use row with [0, s] to simplify the second column. Reorder first.
    M = sympy.Matrix([M.row(0), M.row(2), M.row(1), M.row(3)])
    # R3_new = R3 - s*R2 (in the new order)
    M.row_op(2, lambda val, j: sympy.expand(val - s * M[1, j]))
    # R4_new = R4 + R2 (in the new order)
    M.row_op(3, lambda val, j: sympy.expand(val + M[1, j]))

    # At this point, the matrix is:
    # [ -1,  1 ]
    # [  0,  s ]
    # [  0, -1 ]
    # [  0,  0 ]

    # Step 6: The row [0, -1] is the best pivot for the second column. Reorder and use it.
    M = sympy.Matrix([M.row(0), M.row(2), M.row(1), M.row(3)])
    
    # At this point, the matrix is:
    # [ -1,  1 ]
    # [  0, -1 ]
    # [  0,  s ]
    # [  0,  0 ]

    # Step 7: Perform backward elimination. Use R2 to clear elements in R1 and R3.
    # R1_new = R1 + R2
    M.row_op(0, lambda val, j: sympy.expand(val + M[1, j]))
    # R3_new = R3 + s*R2
    M.row_op(2, lambda val, j: sympy.expand(val + s * M[1, j]))
    
    # The final reduced matrix is:
    # [ -1, 0 ]
    # [  0,-1 ]
    # [  0, 0 ]
    # [  0, 0 ]

    # Step 8: The GCRD is the non-zero part of the reduced matrix.
    G = M[0:2, 0:2]
    
    print("The Greatest Common Right Divisor (GCRD) is the matrix G:")
    g11 = G[0, 0]
    g12 = G[0, 1]
    g21 = G[1, 0]
    g22 = G[1, 1]
    
    print(f"[[{g11}, {g12}],")
    print(f" [{g21}, {g22}]]")

solve()