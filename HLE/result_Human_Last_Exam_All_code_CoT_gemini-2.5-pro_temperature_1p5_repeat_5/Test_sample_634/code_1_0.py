import sympy

def find_gcrd():
    """
    This function calculates the Greatest Common Right Divisor (GCRD)
    of two given polynomial matrices P1 and P2.
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])

    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # Form the stacked matrix M = [P1; P2]
    M = sympy.Matrix.vstack(P1, P2)
    
    # R is a mutable copy of M for row operations
    R = M.copy()

    # Step 1: Use rows from P2 to simplify rows from P1
    # R[0] = R[0] - (s + 1) * R[2]
    # This makes R[0,0] zero
    R[0, :] = R[0, :] - (s + 1) * R[2, :]
    
    # R[1] = R[1] + s * R[2]
    # This simplifies R[1,0] to a constant
    R[1, :] = R[1, :] + s * R[2, :]
    
    # Step 2: Use the new row R[1] (which is [-1, s**2]) as a pivot
    # to eliminate the first column in other rows.
    # Swap R[1] to the top for convenience
    R.row_swap(0, 1)
    
    # R[1] = R[1] + s * R[0]
    R[1, :] = R[1, :] + s * R[0, :]
    
    # R[3] = R[3] - (s + 1) * R[0]
    R[3, :] = R[3, :] - (s + 1) * R[0, :]

    # After these operations, the first column below the first row is zero.
    # The matrix is now:
    # [-1,            s**2]
    # [ 0,            -s**3]  (Original R[0] after ops)
    # [ s,                0]
    # [ 0, 1 - s**3 - s**2]
    # We made a mistake in tracking the rows. Let's operate on the current state.
    # The current matrix state after Step 2 (before re-calculating other rows) is:
    # [      -1,       s**2]
    # [       0,         -s]
    # [       s,          0]
    # [  -s - 1,          1]
    
    # Let's restart the row reduction from a clean state for clarity in the code.
    R = M.copy()
    
    # R1 -> R1 - (s+1)*R3 => [0, -s]
    R[0,:] = R[0,:] - (s+1)*R[2,:]
    # R2 -> R2 + s*R3 => [-1, s^2]
    R[1,:] = R[1,:] + s*R[2,:]
    
    # Now matrix R is:
    # [   0,   -s ]
    # [  -1,  s^2 ]
    # [   s,    0 ]
    # [ -s-1,   1 ]
    
    # Swap R0 and R1 to get a better pivot
    R.row_swap(0,1)
    # Matrix R is:
    # [  -1,  s^2 ]
    # [   0,   -s ]
    # [   s,    0 ]
    # [ -s-1,   1 ]
    
    # Eliminate other elements in column 1
    # R2 -> R2 + s*R0
    R[2,:] = R[2,:] + s*R[0,:]
    # R3 -> R3 - (s+1)*R0
    R[3,:] = R[3,:] - (s+1)*R[0,:]
    # Matrix R is:
    # [  -1,           s**2 ]
    # [   0,             -s ]
    # [   0,            s^3 ]
    # [   0, 1 - s**3 - s**2]

    # Now use R1 = [0, -s] to eliminate column 2 elements
    # R0 -> R0 + s*R1
    R[0,:] = R[0,:] + s*R[1,:]
    # R2 -> R2 + s^2*R1
    R[2,:] = R[2,:] + s**2 * R[1,:]
    # R3 -> R3 - (s^2+s)*R1
    R[3,:] = R[3,:] - (sympy.poly(s**2+s,s)) * R[1,:]
    # Matrix R is:
    # [ -1, 0 ]
    # [  0,-s ]
    # [  0, 0 ]
    # [  0, 1 ]

    # Use R3=[0,1] to eliminate R1
    # R1 -> R1 + s*R3
    R[1,:] = R[1,:] + s*R[3,:]
    # Matrix R is:
    # [ -1, 0 ]
    # [  0, 0 ]
    # [  0, 0 ]
    # [  0, 1 ]

    # Clean up the matrix by swapping rows
    R.row_swap(1,3)
    
    # The final reduced matrix (in echelon form)
    G_unsimplified = R[:2,:]

    # A GCRD is unique up to left multiplication by a unimodular matrix.
    # The matrix [[-1, 0], [0, 1]] is a valid GCRD.
    # For simplicity, we can multiply the first row by -1 to get the identity matrix.
    G = G_unsimplified.copy()
    G[0,:] = -1 * G[0,:]

    # The GCRD is the identity matrix, which means P1 and P2 are right co-prime.
    gcrd_matrix = G
    
    # Print the result in the requested format
    print("The Greatest Common Right Divisor (GCRD) matrix is:")
    print(f"[[{gcrd_matrix[0,0]}, {gcrd_matrix[0,1]}],")
    print(f" [{gcrd_matrix[1,0]}, {gcrd_matrix[1,1]}]]")

find_gcrd()