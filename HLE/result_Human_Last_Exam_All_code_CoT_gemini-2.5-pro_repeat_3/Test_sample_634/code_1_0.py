import sympy

def find_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices
    using elementary row operations on a stacked matrix.
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

    # Form the stacked matrix M by joining P1 and P2 vertically
    M = P1.col_join(P2)

    # The following sequence of elementary row operations reduces M to the form [G; 0]

    # Swap R1 and R3 to get a simpler pivot
    M.row_swap(0, 2)
    # R2 -> R2 + s*R1
    M[1, :] = M[1, :] + s * M[0, :]
    # R3 -> R3 - (s+1)*R1
    M[2, :] = M[2, :] - (s + 1) * M[0, :]

    # Swap R1 and R2 to bring [-1, s**2] to the top, then normalize it to [1, -s**2]
    M.row_swap(0, 1)
    M[0, :] = -1 * M[0, :]

    # Use the new R1 = [1, -s**2] to eliminate the first element in other rows
    # R2 -> R2 - s*R1
    M[1, :] = M[1, :] - s * M[0, :]
    # R4 -> R4 + (s+1)*R1
    M[3, :] = M[3, :] + (s + 1) * M[0, :]
    
    # At this point, M has the form:
    # [1,     -s**2]
    # [0,     s**3]
    # [0,     -s]
    # [0,     1 - s**3 - s**2]

    # We now focus on the second column. Using polynomial division, we know that
    # 1 = (1 - s**3 - s**2) - (s**2 + s)*(-s).
    # This translates to the row operation R4_new = R4 - (s**2 + s)*R3, which will produce [0, 1].
    M[3, :] = M[3, :] - (s**2 + s) * M[2, :]
    
    # Swap R2 and R4 to use [0, 1] as a pivot
    M.row_swap(1, 3)

    # Finally, use the pivot row R2 = [0, 1] to zero out the second column in other rows
    # R1 -> R1 + s**2 * R2
    M[0, :] = M[0, :] + s**2 * M[1, :]
    # R3 -> R3 - s**3 * R2 (M[2,:] is now [0, s**3])
    M[2, :] = M[2, :] - s**3 * M[1, :]
    # R4 -> R4 + s * R2 (M[3,:] is now [0, -s])
    M[3, :] = M[3, :] + s * M[1, :]

    # The top 2x2 matrix is the GCRD
    GCRD = M[:2, :]

    print("The greatest common right divisor is the matrix G(s) =")
    print(f"[[{GCRD[0,0]}, {GCRD[0,1]}];")
    print(f" [{GCRD[1,0]}, {GCRD[1,1]}]]")


if __name__ == "__main__":
    find_gcrd()