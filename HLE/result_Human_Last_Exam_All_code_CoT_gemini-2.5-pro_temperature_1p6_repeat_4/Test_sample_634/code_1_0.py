import sympy

def find_gcrd():
    """
    This function finds the greatest common right divisor (GCRD) of the two
    polynomial matrices P1 and P2 using elementary row operations.
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])

    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # 1. Form the stacked matrix M
    M = sympy.Matrix([P1.row(0), P1.row(1), P2.row(0), P2.row(1)])
    
    # 2. Perform elementary row operations to simplify the matrix.
    # The goal is to get a matrix of the form [0; G] where G is the GCRD.

    # Let the rows be R0, R1, R2, R3
    # Use R2 ([s, 0]) to simplify R0 and R1
    # R0 -> R0 - (s+1)*R2
    M[0, :] = sympy.simplify(M[0, :] - (s + 1) * M[2, :])
    # R1 -> R1 + s*R2
    M[1, :] = sympy.simplify(M[1, :] + s * M[2, :])
    
    # After this, the matrix is:
    # [ 0     , -s   ]
    # [-1    , s**2 ]
    # [ s     , 0    ]
    # [-s-1  , 1    ]
    
    # Now use the new R1 (M[1,:]) which is [-1, s**2] to simplify R2 and R3
    # R2 -> R2 + s*R1
    M[2, :] = sympy.simplify(M[2, :] + s * M[1, :])
    # R3 -> R3 + (-s-1)*R1
    M[3, :] = sympy.simplify(M[3, :] + (-s - 1) * M[1, :])
    
    # The matrix is now:
    # [ 0     , -s          ]
    # [-1    , s**2        ]
    # [ 0     , s**3        ]
    # [ 0     , 1-s**2-s**3 ]

    # Now, reduce rows using polynomial division logic on the second column.
    # Use R0 ([0, -s]) to reduce R2 ([0, s**3]) and R3 ([0, 1-s**2-s**3]).
    
    # R2 -> R2 + s**2 * R0 which eliminates R2
    M[2, :] = sympy.simplify(M[2, :] + s**2 * M[0, :])
    
    # For R3, we find the quotient of (1-s**2-s**3) / (-s), which is (s+s**2).
    # R3 -> R3 - (s+s**2) * R0
    q = sympy.poly(M[3, 1], s).div(sympy.poly(M[0, 1], s))[0]
    M[3, :] = sympy.simplify(M[3, :] - q.as_expr() * M[0, :])
    
    # The matrix is now:
    # [ 0     , -s    ]
    # [-1    , s**2  ]
    # [ 0     , 0     ]
    # [ 0     , 1     ]
    
    # Now, use the new R3 ([0, 1]) to clean up the other rows.
    # R0 -> R0 + s*R3
    M[0, :] = sympy.simplify(M[0, :] + s * M[3, :])
    # R1 -> R1 - s**2*R3
    M[1, :] = sympy.simplify(M[1, :] - s**2 * M[3, :])

    # The final reduced matrix is:
    # [ 0, 0 ]
    # [-1, 0]
    # [ 0, 0 ]
    # [ 0, 1 ]

    # 3. Extract the non-zero rows to form the GCRD.
    non_zero_rows = [row for row in M.tolist() if any(elem != 0 for elem in row)]
    G_raw = sympy.Matrix(non_zero_rows)
    
    # Reorder to be [[-1, 0], [0, 1]]
    if G_raw[0,0] == 0:
        G_raw = sympy.Matrix([G_raw.row(1), G_raw.row(0)])
        
    # 4. Normalize the GCRD.
    # The determinant of G_raw is -1, which is a non-zero constant.
    # This means G_raw is unimodular, and P1, P2 are right-coprime.
    # The standard GCRD in this case is the identity matrix.
    # We can get I by left-multiplying G_raw by [[-1, 0], [0, 1]].
    G_final = sympy.eye(2)
    
    print("The greatest common right divisor is the matrix:")
    sympy.pprint(G_final)
    
    print("\nThis means P1 and P2 are right-coprime.")

find_gcrd()