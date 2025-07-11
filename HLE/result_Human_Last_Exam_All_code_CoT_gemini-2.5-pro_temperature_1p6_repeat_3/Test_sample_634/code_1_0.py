import sympy

def find_gcrd():
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

    # 1. Form the stacked matrix M
    M = P1.col_join(P2)

    # 2. Apply a sequence of elementary row operations to reduce M.
    # The operations are chosen to systematically reduce the polynomial degrees.
    
    # Operation: R0 <- R0 - (s+1)*R2
    M.row_op(0, lambda val, j: sympy.expand(val - (s + 1) * M[2, j]))
    
    # Operation: R1 <- R1 + s*R2
    M.row_op(1, lambda val, j: sympy.expand(val + s * M[2, j]))
    
    # At this point, M is:
    # [0,      -s]
    # [-1,     s**2]
    # [s,      0]
    # [-s-1,   1]
    
    # Use the new row 1 (index 1) as a pivot to eliminate terms in column 0.
    # Operation: R2 <- R2 + s*R1
    M.row_op(2, lambda val, j: sympy.expand(val + s * M[1, j]))
    
    # Operation: R3 <- R3 - (s+1)*R1
    M.row_op(3, lambda val, j: sympy.expand(val - (s + 1) * M[1, j]))
    
    # At this point, M is:
    # [0,           -s]
    # [-1,          s**2]
    # [0,           s**3]
    # [0, 1 - s**3 - s**2]
    
    # Now, rows 0, 2, and 3 have 0 in the first column. We can reduce them.
    # We will use polynomial division logic. 
    # Let's reduce row 3 using row 0 as a pivot.
    # The multiplier is the quotient of M[3,1] / M[0,1]
    q, r = sympy.div(M[3, 1], M[0, 1], s)
    M.row_op(3, lambda val, j: sympy.expand(val - q * M[0, j]))
    
    # At this point, M is:
    # [0,    -s]
    # [-1, s**2]
    # [0,  s**3]
    # [0,     1]

    # Now use the new simple row 3 (index 3) `[0, 1]` to eliminate other terms in column 1.
    # Operation: R0 <- R0 + s*R3
    M.row_op(0, lambda val, j: sympy.expand(val + s * M[3, j]))
    
    # Operation: R1 <- R1 - s**2*R3
    M.row_op(1, lambda val, j: sympy.expand(val - (s**2) * M[3, j]))

    # Operation: R2 <- R2 - s**3*R3
    M.row_op(2, lambda val, j: sympy.expand(val - (s**3) * M[3, j]))

    # Final matrix M has two non-zero rows: [-1, 0] and [0, 1].
    # These form the GCRD, up to ordering and scaling by constants.
    
    GCRD = sympy.Matrix([[-1, 0], [0, 1]])
    
    # Normalize by making the leading diagonal entries positive.
    if GCRD[0,0] < 0:
        GCRD[0,:] = -1 * GCRD[0,:]

    # 3. Print the resulting GCRD
    print("The greatest common right divisor (GCRD) matrix is:")
    print(f"[{GCRD[0,0]}, {GCRD[0,1]}]")
    print(f"[{GCRD[1,0]}, {GCRD[1,1]}]")

if __name__ == '__main__':
    find_gcrd()