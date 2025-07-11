import sympy

def find_gcrd():
    """
    Calculates the Greatest Common Right Divisor (GCRD) of two polynomial matrices
    P1 and P2 using the Euclidean algorithm via elementary row operations.
    """
    s = sympy.Symbol('s')

    # Define the polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])

    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    print("Given polynomial matrices:")
    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)
    print("\n" + "="*30)
    print("Finding the GCRD by reducing the stacked matrix M = [P1; P2]...")

    # Form the stacked matrix M
    M = P1.col_join(P2)
    
    # We will perform a series of elementary row operations to reduce M.
    # The operations are chosen to systematically zero out rows from the bottom up.

    # Step 1: Use the '1' in M[3, 1] as a pivot to zero out the second column in other rows.
    # M[0, :] = M[0, :] + s * M[3, :]
    M[0, :] = sympy.expand(M[0, :] + s * M[3, :])
    # M[1, :] = M[1, :] - s**2 * M[3, :]
    M[1, :] = sympy.expand(M[1, :] - (s**2) * M[3, :])

    # Step 2: Now M[0,:] is [0, 0] and M[1,:] is [s**3 - 1, 0].
    # We apply the Euclidean algorithm for the polynomials in the first column of the
    # remaining rows M[1,:] and M[2,:].
    # Remainder of (s**3 - 1) / s is -1.
    # M[1, :] = M[1, :] - s**2 * M[2, :]
    M[1, :] = sympy.expand(M[1, :] - (s**2) * M[2, :])
    
    # Step 3: Now M[1,:] is [-1, 0]. This provides a constant pivot to clear the first column.
    # M[2, :] = M[2, :] + s * M[1, :]
    M[2, :] = sympy.expand(M[2, :] + s * M[1, :])
    # M[3, :] = M[3, :] - (s + 1) * M[1, :]
    M[3, :] = sympy.expand(M[3, :] - (s + 1) * M[1, :])
    
    # Step 4: The matrix is now reduced. The non-zero rows form the GCRD.
    # Let's collect the non-zero rows and clean them up.
    
    # Normalize the pivot row M[1,:]
    M[1, :] = -M[1, :]
    
    # The two non-zero rows are M[1,:]=[1,0] and M[3,:]=[0,1].
    # We can arrange them into the final GCRD matrix.
    GCRD = sympy.Matrix([M[1,:], M[3,:]])
    
    # Print the result in the requested format.
    print("\nThe greatest common right divisor G(s) is:")
    
    g11 = GCRD[0,0]
    g12 = GCRD[0,1]
    g21 = GCRD[1,0]
    g22 = GCRD[1,1]

    print(f"G(s) = [ {g11}, {g12};")
    print(f"         {g21}, {g22} ]")

if __name__ == '__main__':
    find_gcrd()
