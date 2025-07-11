import sympy

def find_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two 
    polynomial matrices P1 and P2 using elementary row operations.
    """
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

    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)
    print("\n" + "="*30 + "\n")

    # Step 1: Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("Step 1: Form the stacked matrix M = [P1; P2]")
    sympy.pprint(M)
    print("-" * 30)

    # Step 2: Perform row reduction to find the GCRD
    # Create a copy to modify
    M_reduced = M.copy()

    # R1 <-> R3 to get a simpler pivot
    M_reduced.row_swap(0, 2)
    # R2_new = R2 + s*R1
    M_reduced[1, :] += s * M_reduced.row(0)
    # R3_new = R3 - (s+1)*R1
    M_reduced[2, :] -= (s + 1) * M_reduced.row(0)
    # R4_new = R4 + R1
    M_reduced[3, :] += M_reduced.row(0)
    
    # Now M_reduced is [[s, 0], [-1, s^2], [0, -s], [-1, 1]]
    # R1 <-> R4 to get a constant pivot
    M_reduced.row_swap(0, 3)
    
    # Now M_reduced is [[-1, 1], [-1, s^2], [0, -s], [s, 0]]
    # R2_new = R2 - R1
    M_reduced[1, :] -= M_reduced.row(0)
    # R4_new = R4 + s*R1
    M_reduced[3, :] += s * M_reduced.row(0)

    # Now M_reduced is [[-1, 1], [0, s^2-1], [0, -s], [0, s]]
    # R3_new = R3 + R4
    M_reduced[2, :] += M_reduced.row(3)

    # Now M_reduced is [[-1, 1], [0, s^2-1], [0, 0], [0, s]]
    # R2_new = R2 - s*R4
    M_reduced[1, :] -= s * M_reduced.row(3)

    # Now M_reduced is [[-1, 1], [0, -1], [0, 0], [0, s]]
    # R1_new = R1 + R2
    M_reduced[0, :] += M_reduced.row(1)
    # R4_new = R4 + s*R2
    M_reduced[3, :] += s * M_reduced.row(1)

    # Now M_reduced is [[-1, 0], [0, -1], [0, 0], [0, 0]]
    print("Step 2: After performing elementary row operations, the matrix becomes:")
    sympy.pprint(M_reduced)
    print("-" * 30)

    # Step 3: The non-zero rows form the GCRD. Normalize them.
    GCRD = sympy.Matrix([
        [-1 * M_reduced.row(0)],
        [-1 * M_reduced.row(1)]
    ]).applyfunc(sympy.simplify)

    print("Step 3: The non-zero rows, normalized, form the GCRD.")
    print("\nThe Greatest Common Right Divisor is:")
    sympy.pprint(GCRD)

    # Print the final answer in the requested format
    print("\nFinal Answer in equation form:")
    print(f"GCRD = [[{GCRD[0,0]}, {GCRD[0,1]}], [{GCRD[1,0]}, {GCRD[1,1]}]]")

if __name__ == '__main__':
    find_gcrd()