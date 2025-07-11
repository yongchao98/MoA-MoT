import sympy as sp

def find_gcrd():
    """
    Finds the Greatest Common Right Divisor (GCRD) of two polynomial matrices P1 and P2.
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

    print("Given polynomial matrices:")
    print("P1(s) =")
    sp.pprint(P1)
    print("\nP2(s) =")
    sp.pprint(P2)

    # Step 1: Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)
    print("\nStep 1: Form the stacked matrix M = [P1; P2]")
    sp.pprint(M)

    # Step 2: Reorder rows to have lower-degree polynomials at the top.
    # Swap R1 <-> R3 and R2 <-> R4 for easier computation.
    M.row_swap(0, 2)
    M.row_swap(1, 3)
    print("\nStep 2: Reorder rows to have lower degrees on top")
    sp.pprint(M)

    # Step 3: Perform elementary row operations to introduce zeros.
    print("\nStep 3: Perform row reduction to find the GCRD.")
    # The '1' in the new R2 (M[1,:]) is useful for elimination.

    # R3 -> R3 + s*R2
    # [s**2+s, -s] + s*[-s-1, 1] = [s**2+s - s**2-s, -s+s] = [0, 0]
    M[2, :] = sp.simplify(M[2, :] + s * M[1, :])
    print("\nAfter R3 -> R3 + s*R2:")
    sp.pprint(M)

    # R4 -> R4 - s**2 * R2
    # [-s**2-1, s**2] - s**2*[-s-1, 1] = [-s**2-1 - (-s**3-s**2), s**2-s**2] = [s**3-1, 0]
    M[3, :] = sp.simplify(M[3, :] - s**2 * M[1, :])
    print("\nAfter R4 -> R4 - s^2*R2:")
    sp.pprint(M)

    # R4 -> R4 - s**2 * R1
    # [s**3-1, 0] - s**2*[s, 0] = [s**3-1-s**3, 0] = [-1, 0]
    M[3, :] = sp.simplify(M[3, :] - s**2 * M[0, :])
    print("\nAfter R4 -> R4 - s^2*R1:")
    sp.pprint(M)

    # Swap R1 and R4 to get a constant in the top-left position
    M.row_swap(0, 3)
    print("\nAfter swapping R1 and R4:")
    sp.pprint(M)

    # R4 -> R4 + s*R1
    # [s, 0] + s*[-1, 0] = [0, 0]
    M[3, :] = sp.simplify(M[3, :] + s * M[0, :])
    print("\nAfter R4 -> R4 + s*R1:")
    sp.pprint(M)

    # The top 2x2 matrix is a GCRD
    G = M[:2, :]
    print("\nStep 4: A GCRD is the non-zero part of the matrix:")
    sp.pprint(G)

    # Step 5: Normalize the GCRD to its simplest form (Hermite Normal Form)
    # R1 -> -1 * R1 to make the pivot 1.
    G[0, :] = -1 * G[0, :]

    # R2 -> R2 + (s+1) * R1 to make the matrix upper-triangular.
    G[1, :] = sp.simplify(G[1, :] + (s+1) * G[0, :])
    print("\nStep 5: Normalize the GCRD to its simplest form (Identity Matrix):")
    sp.pprint(G)

    # Final result output
    print("\n" + "="*50)
    print("The final Greatest Common Right Divisor (GCRD) is the identity matrix:")
    print("G = [[1, 0],")
    print("     [0, 1]]")
    print("\nThe elements of the final GCRD matrix G are:")
    print(f"G[0, 0] = {G[0,0]}")
    print(f"G[0, 1] = {G[0,1]}")
    print(f"G[1, 0] = {G[1,0]}")
    print(f"G[1, 1] = {G[1,1]}")
    print("="*50)

if __name__ == '__main__':
    find_gcrd()
