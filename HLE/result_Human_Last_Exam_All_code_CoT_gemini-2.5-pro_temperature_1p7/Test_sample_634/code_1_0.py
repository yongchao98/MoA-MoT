import sympy

def find_gcrd():
    """
    This function calculates the Greatest Common Right Divisor (GCRD) of two
    given polynomial matrices P1 and P2 by performing elementary row operations
    on their stacked form.
    """
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

    # 1. Stack P1 and P2 to form matrix M
    M = sympy.Matrix.vstack(P1, P2)
    
    print("Initial stacked matrix M = [P1; P2]:")
    sympy.pprint(M)
    print("-" * 30)

    # --- Start of Row Reduction ---
    # The process follows the standard Gaussian elimination algorithm adapted for polynomial matrices.
    
    # Swap R1 and R3 to get a simpler pivot (degree 1) at the top.
    M.row_swap(0, 2)
    
    # Use the pivot row R1 ([s, 0]) to simplify other rows.
    # R2 -> R2 + s*R1
    M[1,:] = sympy.poly(M[1,:] + s * M[0,:], s).as_expr()
    # R3 -> R3 - (s+1)*R1
    M[2,:] = sympy.poly(M[2,:] - (s+1) * M[0,:], s).as_expr()
    # R4 -> R4 + R1
    M[3,:] = sympy.poly(M[3,:] + M[0,:], s).as_expr()

    # The matrix is now:
    # [ s,   0    ]
    # [ -1,  s^2  ]
    # [ 0,  -s   ]
    # [ -1,  1    ]
    
    # Use R4 ([-1, 1]) as the next pivot. Swap it up to R2.
    M.row_swap(1, 3)
    # Make the pivot element 1 by multiplying the row by -1.
    M[1,:] = -1 * M[1,:]
    
    # Use the new pivot R2 ([1, -1]) to eliminate other entries in the first column.
    # R1 -> R1 - s*R2
    M[0,:] = sympy.poly(M[0,:] - s * M[1,:], s).as_expr()
    # R4 -> R4 + R2
    M[3,:] = sympy.poly(M[3,:] + M[1,:], s).as_expr()
    
    # The matrix is now:
    # [ 0,   s      ]
    # [ 1,  -1      ]
    # [ 0,  -s      ]
    # [ 0,  s^2 - 1 ]
    
    # Reorder rows for a standard echelon form.
    M.row_swap(0, 1)
    
    # Use R2 ([0, s]) to simplify rows below it.
    # R3 -> R3 + R2
    M[2,:] = sympy.poly(M[2,:] + M[1,:], s).as_expr()
    # R4 -> R4 - s*R2
    M[3,:] = sympy.poly(M[3,:] - s * M[1,:], s).as_expr()
    
    # The matrix is now:
    # [ 1,  -1   ]
    # [ 0,   s   ]
    # [ 0,   0   ]
    # [ 0,  -1   ]
    
    # Use R4 ([0, -1]) as the next pivot.
    M.row_swap(1, 3)
    # Normalize the pivot to 1.
    M[1,:] = -1 * M[1,:]
    # Swap zero rows to the bottom.
    M.row_swap(2, 3)
    
    # Use the new pivot R2 ([0, 1]) to simplify R3 ([0, s]).
    # R3 -> R3 - s*R2
    M[2,:] = sympy.poly(M[2,:] - s * M[1,:], s).as_expr()

    # --- End of Row Reduction ---
    
    # The final reduced matrix is:
    # [ 1,  -1   ]
    # [ 0,   1   ]
    # [ 0,   0   ]
    # [ 0,   0   ]
    
    # The non-zero rows form the GCRD.
    non_zero_rows = []
    for r in range(M.rows):
        if not M.row(r).is_zero:
            non_zero_rows.append(M.row(r))
            
    GCRD = sympy.Matrix(non_zero_rows)
    
    # The final step is to put the GCRD into its reduced row-echelon form.
    # Use the pivot in R2 to eliminate the entry above it in R1.
    # R1 -> R1 + R2
    GCRD[0,:] = GCRD[0,:] + GCRD[1,:]

    print("The final GCRD in reduced row-echelon form is:")
    print("GCRD = ")
    sympy.pprint(GCRD)
    
    # Printing each number in the final matrix explicitly.
    final_matrix_list = GCRD.tolist()
    print("\nIn list format, the GCRD is:")
    print(final_matrix_list)


if __name__ == '__main__':
    find_gcrd()
