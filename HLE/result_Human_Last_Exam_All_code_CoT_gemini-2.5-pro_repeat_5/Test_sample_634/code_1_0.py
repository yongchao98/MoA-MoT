import sympy

def find_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial
    matrices, P1 and P2, using the Euclidean algorithm.
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

    print("Given polynomial matrices:")
    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)

    # Form the stacked matrix M = [P1; P2] for the Euclidean algorithm.
    # The rows of M are indexed 0, 1, 2, 3.
    M = P1.vstack(P2)
    print("\nStarting the reduction process on the composite matrix M = [P1; P2]:")
    
    # --- Step 1: Create a simple pivot row ---
    # Operation: R3 -> R3 + R4 (M[2,:] -> M[2,:] + M[3,:])
    # This creates the pivot row [-1, 1]
    M[2, :] = M[2, :] + M[3, :]

    # --- Step 2: First round of eliminations ---
    # Use the new R3 = [-1, 1] as a pivot to eliminate the first element in other rows.
    M[0, :] = M[0, :] + (s**2 + s) * M.row(2)
    M[1, :] = M[1, :] - (s**2 + 1) * M.row(2)
    M[3, :] = M[3, :] - (s + 1) * M.row(2)
    
    # Simplify the matrix after the operations
    M = sympy.simplify(M)
    
    # --- Step 3: Second round of eliminations ---
    # At this point, M has rows [0, s**2], [0, -1], [-1, 1], [0, -s].
    # We will use R2 (now [0, -1]) as the next pivot to eliminate the second element in other rows.
    M[0, :] = M[0, :] + (s**2) * M.row(1)
    M[2, :] = M[2, :] + M.row(1)
    M[3, :] = M[3, :] - s * M.row(1)

    # Simplify the final matrix
    M = sympy.simplify(M)
    
    # --- Step 4: Extract the GCRD ---
    # The non-zero rows of the final matrix form the GCRD.
    gcrd_rows = [row for row in M.tolist() if not all(x == 0 for x in row)]
    GCRD = sympy.Matrix(gcrd_rows)

    print("\nThe final reduced matrix (with zero rows removed) is the GCRD:")
    sympy.pprint(GCRD)
    
    print("\nThe elements of the GCRD matrix G(s) are:")
    print(f"G = [ {GCRD[0,0]}, {GCRD[0,1]} ]")
    print(f"    [ {GCRD[1,0]}, {GCRD[1,1]} ]")

if __name__ == "__main__":
    find_gcrd()