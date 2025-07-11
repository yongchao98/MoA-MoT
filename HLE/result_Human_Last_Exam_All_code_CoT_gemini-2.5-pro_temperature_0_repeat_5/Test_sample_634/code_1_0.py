import sympy as sp

def solve_gcrd():
    """
    This function calculates the Greatest Common Right Divisor (GCRD)
    of two polynomial matrices P1 and P2.
    """
    # Define the symbolic variable 's'
    s = sp.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sp.Matrix([[s**2 + s, -s],
                    [-s**2 - 1, s**2]])

    P2 = sp.Matrix([[s, 0],
                    [-s - 1, 1]])

    # Step 1: Form the augmented matrix M by stacking P1 and P2
    M = P1.col_join(P2)
    print("Initial augmented matrix M = [P1; P2]:")
    sp.pprint(M)
    print("-" * 30)

    # Step 2: Perform elementary row operations to reduce M.
    # Let the rows be R0, R1, R2, R3.
    
    # Operation: R3 -> R3 + R2
    # This creates a row with constant entries, which is an excellent pivot.
    M[3, :] = M[3, :] + M[2, :]
    print("After R3 -> R3 + R2:")
    sp.pprint(M)
    print("-" * 30)
    
    # Now, use the new R3 = [-1, 1] to eliminate the first column in other rows.
    # Operation: R0 -> R0 + (s**2 + s) * R3
    M[0, :] = M[0, :] + (s**2 + s) * M[3, :]
    # Operation: R1 -> R1 - (s**2 + 1) * R3
    M[1, :] = M[1, :] - (s**2 + 1) * M[3, :]
    # Operation: R2 -> R2 + s * R3
    M[2, :] = M[2, :] + s * M[3, :]
    
    print("After using R3 to eliminate other entries in the first column:")
    sp.pprint(M)
    print("-" * 30)

    # The matrix is now:
    # [  0,   s**2]
    # [  0,     -1]
    # [  0,      s]
    # [ -1,      1]
    # Let's use R1 = [0, -1] as the pivot for the second column.
    
    # Operation: R0 -> R0 + s**2 * R1
    M[0, :] = M[0, :] + (s**2) * M[1, :]
    # Operation: R2 -> R2 + s * R1
    M[2, :] = M[2, :] + s * M[1, :]
    # Operation: R3 -> R3 + 1 * R1
    M[3, :] = M[3, :] + M[1, :]

    print("After using R1 to eliminate other entries in the second column:")
    sp.pprint(M)
    print("-" * 30)

    # The matrix is now:
    # [  0, 0]
    # [  0,-1]
    # [  0, 0]
    # [ -1, 0]
    # The non-zero rows are [-1, 0] and [0, -1].
    
    # Step 3: Extract the non-zero rows to form the GCRD.
    GCRD = sp.Matrix([M[3, :], M[1, :]])

    print("The final non-zero rows form the GCRD matrix.")
    print("Greatest Common Right Divisor (GCRD):")
    
    # Print the final matrix elements as requested
    row1 = GCRD.row(0)
    row2 = GCRD.row(1)
    print(f"[{int(row1[0])}, {int(row1[1])}]")
    print(f"[{int(row2[0])}, {int(row2[1])}]")

if __name__ == '__main__':
    solve_gcrd()