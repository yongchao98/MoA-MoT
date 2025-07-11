import sympy as sp

def find_gcrd():
    """
    This function calculates the greatest common right divisor (GCRD) of
    two given polynomial matrices P1 and P2.
    """
    # Define the symbolic variable 's'
    s = sp.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sp.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])
    P2 = sp.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    # 1. Form the stacked matrix M
    M = P1.col_join(P2)

    # 2. Perform elementary row operations to reduce M.
    # We will follow a systematic procedure similar to Gaussian elimination.

    # Make the pivot in the first column simpler by swapping R1 and R3
    M[0, :], M[2, :] = M[2, :], M[0, :]

    # Use the new R1 = [s, 0] to simplify other rows
    M[1, :] = sp.expand(M.row(1) + s * M.row(0))
    M[2, :] = sp.expand(M.row(2) - (s + 1) * M.row(0))
    M[3, :] = sp.expand(M.row(3) + M.row(0))
    # M is now:
    # [ s,    0   ]
    # [ -1,   s**2 ]
    # [ 0,    -s  ]
    # [ -1,   1   ]

    # Make the pivot -1 by swapping R1 and R2
    M[0, :], M[1, :] = M[1, :], M[0, :]

    # Use the new R1 = [-1, s**2] to eliminate other entries in the first column
    M[1, :] = sp.expand(M.row(1) + s * M.row(0))
    M[3, :] = sp.expand(M.row(3) - M.row(0))
    # M is now:
    # [ -1,   s**2   ]
    # [ 0,    s**3   ]
    # [ 0,    -s    ]
    # [ 0,    1-s**2 ]

    # Now focus on the second column for rows 2, 3, 4.
    # Use R3 = [0, -s] as pivot to simplify R2 and R4.
    M[1, :] = sp.expand(M.row(1) + s**2 * M.row(2))
    M[3, :] = sp.expand(M.row(3) - s * M.row(2))
    # M is now:
    # [ -1, s**2 ]
    # [ 0,    0  ]
    # [ 0,   -s  ]
    # [ 0,    1  ]

    # Use R4 = [0, 1] as the pivot for the second column
    # Swap R2 and R4 to bring the pivot up
    M[1, :], M[3, :] = M[3, :], M[1, :]

    # Use the new R2 = [0, 1] to eliminate other entries in the second column
    M[0, :] = sp.expand(M.row(0) - s**2 * M.row(1))
    M[2, :] = sp.expand(M.row(2) + s * M.row(1))
    # M is now:
    # [ -1, 0 ]
    # [ 0,  1 ]
    # [ 0,  0 ]
    # [ 0,  0 ]

    # 3. Extract the non-zero rows which form the GCRD matrix G.
    G = M[0:2, :]

    # 4. Print the resulting GCRD matrix.
    print("The greatest common right divisor G(s) is:")
    print(f"G(s) = [{G[0,0]}, {G[0,1]};")
    print(f"       {G[1,0]},  {G[1,1]}]")


if __name__ == "__main__":
    find_gcrd()
