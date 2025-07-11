import sympy

def find_gcrd():
    """
    This function calculates the greatest common right divisor (GCRD) of
    two given polynomial matrices P1 and P2.
    """
    # Define the symbolic variable 's'
    s = sympy.symbols('s')

    # Define the polynomial matrices P1 and P2 as given in the problem
    P1 = sympy.Matrix([
        [s**2 + s, -s],
        [-s**2 - 1, s**2]
    ])

    P2 = sympy.Matrix([
        [s, 0],
        [-s - 1, 1]
    ])

    print("The given matrices are:")
    print("P1(s) =")
    sympy.pprint(P1)
    print("\nP2(s) =")
    sympy.pprint(P2)
    
    # To find the GCRD, we form a stacked matrix M = [P1; P2]
    # and find its row-echelon form.
    M = P1.col_join(P2)

    # The .rref() method computes the row-reduced echelon form.
    # The non-zero rows of the resulting matrix form the GCRD.
    M_rref = M.rref(pivots=False)

    # Extract the non-zero rows which form the GCRD matrix G.
    first_zero_row = -1
    for i in range(M_rref.rows):
        if M_rref.row(i).is_zero:
            first_zero_row = i
            break
    
    if first_zero_row == -1: # Case where there are no zero rows
        G = M_rref
    else:
        G = M_rref[:first_zero_row, :]

    print("\nThe greatest common right divisor G(s) is:")
    
    # Print the resulting GCRD matrix, showing each element of the "equation"
    for i in range(G.rows):
        row_elements = []
        for j in range(G.cols):
            row_elements.append(str(G[i, j]))
        print(f"[{'   '.join(row_elements)}]")

# Execute the function
find_gcrd()