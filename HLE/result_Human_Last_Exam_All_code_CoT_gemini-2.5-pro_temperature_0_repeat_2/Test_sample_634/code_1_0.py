import sympy

def find_greatest_common_right_divisor():
    """
    This function calculates the greatest common right divisor (GCRD) of two
    given polynomial matrices P1 and P2 using the Euclidean algorithm.
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

    # Form the stacked matrix M = [P1; P2]
    M = P1.col_join(P2)

    # --- Begin Euclidean Algorithm via Elementary Row Operations ---

    # To simplify the process, we want a row with a constant pivot.
    # We can achieve this by adding row 4 to row 3.
    # M[2, :] corresponds to the 3rd row (0-indexed).
    # Operation: R3 -> R3 + R4
    M[2, :] = M[2, :] + M[3, :]

    # Now, row 3 is [-1, 1]. This is an excellent pivot. Let's move it to the top.
    # Operation: Swap R1 and R3
    M.row_swap(0, 2)

    # Now, use the new pivot in R1 (which is [-1, 1]) to create zeros
    # in the first column of the other rows.
    # The pivot element is M[0, 0] = -1.
    # Operation R2 -> R2 - (M[1,0] / M[0,0]) * R1
    M[1, :] = M[1, :] - (M[1, 0] / M[0, 0]) * M[0, :]
    # Operation R3 -> R3 - (M[2,0] / M[0,0]) * R1
    M[2, :] = M[2, :] - (M[2, 0] / M[0, 0]) * M[0, :]
    # Operation R4 -> R4 - (M[3,0] / M[0,0]) * R1
    M[3, :] = M[3, :] - (M[3, 0] / M[0, 0]) * M[0, :]
    
    # After these operations, the matrix has zeros in the first column below the first row.
    # The second row is now [0, -1], which is our next pivot.
    # Use the pivot in R2 (M[1,1] = -1) to create zeros in the second column.
    # Operation R1 -> R1 - (M[0,1] / M[1,1]) * R2
    M[0, :] = M[0, :] - (M[0, 1] / M[1, 1]) * M[1, :]
    # Operation R3 -> R3 - (M[2,1] / M[1,1]) * R2
    M[2, :] = M[2, :] - (M[2, 1] / M[1, 1]) * M[1, :]
    # Operation R4 -> R4 - (M[3,1] / M[1,1]) * R2
    M[3, :] = M[3, :] - (M[3, 1] / M[1, 1]) * M[1, :]

    # The matrix is now reduced. The GCRD is the top 2x2 non-zero block.
    # We need to simplify the expressions in the matrix.
    M = sympy.simplify(M)
    GCRD = M[:2, :]

    # The GCRD is unique up to left multiplication by a unimodular matrix.
    # The result is [[-1, 0], [0, -1]]. We can multiply by [[-1, 0], [0, -1]]
    # (which is unimodular) to get the identity matrix, a simpler form.
    normalizer = sympy.Matrix([[-1, 0], [0, -1]])
    GCRD_normalized = normalizer * GCRD

    # Print the final result
    final_gcrd = GCRD_normalized
    print("The greatest common right divisor is:")
    print(f"[[{final_gcrd[0,0]}, {final_gcrd[0,1]}],")
    print(f" [{final_gcrd[1,0]}, {final_gcrd[1,1]}]]")

if __name__ == '__main__':
    find_greatest_common_right_divisor()