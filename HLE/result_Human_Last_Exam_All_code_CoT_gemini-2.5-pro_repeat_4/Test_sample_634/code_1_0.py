import sympy as sp

def solve():
    """
    Calculates the greatest common right divisor (GCRD) of two polynomial matrices.
    """
    # Define the symbolic variable 's'
    s = sp.Symbol('s')

    # Define the input polynomial matrices P1 and P2
    P1 = sp.Matrix([[s**2 + s, -s],
                    [-s**2 - 1, s**2]])

    P2 = sp.Matrix([[s, 0],
                    [-s - 1, 1]])

    # 1. Form the augmented matrix M by stacking P1 and P2
    M = P1.vstack(P2)

    # 2. Perform elementary row operations to reduce M
    # The process is similar to Gaussian elimination but for polynomials.

    # Step A: Use R2 ([s, 0]) to simplify R0 and R1
    # R0 -> R0 - (s+1)*R2
    M[0, :] = M[0, :] - (s + 1) * M[2, :]
    # R1 -> R1 + s*R2
    M[1, :] = M[1, :] + s * M[2, :]

    # Step B: Simplify the matrix further.
    # R2 -> R2 + R3
    M[2, :] = M[2, :] + M[3, :]
    
    # At this point, M has a simple row M[2,:] = [-1, 1]. We use it as a pivot.
    
    # Step C: Use the new R2 ([-1, 1]) to simplify R1 and R3
    # R1 -> R1 - R2
    M[1, :] = M[1, :] - M[2, :]
    # R3 -> R3 - (s+1)*R2. The factor is chosen to make M[3,0] zero.
    M[3, :] = M[3, :] - (s + 1) * M[2, :]

    # Now the matrix has two zero columns in the first element,
    # and is ready for further reduction. Let's reorder for clarity.
    M.row_swap(0, 2)
    M.row_swap(1, 3)
    
    # Step D: Reduce the second column
    # R2 -> R2 + s*R1
    M[2, :] = M[2, :] + s * M[1, :]
    
    # Step E: Final reduction to get to echelon form
    # R1 -> R1 - s*R2
    M[1, :] = M[1, :] - s * M[2, :]

    # 3. Extract the non-zero rows which form the GCRD
    # The matrix is now in row echelon form, with the bottom two rows being zero.
    G = M[:2, :]
    
    # Normalize the GCRD to its simplest form (an identity matrix in this case)
    # Multiply the first row by -1
    if G[0, 0] == -1:
        G[0, :] = -G[0, :]
    # Multiply the second row by -1
    if G[1, 1] == -1:
        G[1, :] = -G[1, :]

    # Print the final result
    print("The greatest common right divisor is the matrix G(s):")
    # Using a loop to print each element as requested
    for i in range(G.rows):
        row_str = []
        for j in range(G.cols):
            row_str.append(str(G[i, j]))
        print(f"[{', '.join(row_str)}]")

solve()