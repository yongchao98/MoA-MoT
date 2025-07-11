import sympy

def solve_factorization():
    """
    Calculates and verifies a left coprime factorization of a given transfer function matrix H(s).
    The factorization is of the form H(s) = D(s)^-1 * N(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    print("The given transfer function is:")
    print("H(s) =")
    sympy.pretty_print(H)
    print("\n" + "="*40 + "\n")

    # Step 1: Propose a left coprime factorization D(s) and N(s).
    # This result is derived analytically by ensuring the degree of det(D(s))
    # matches the McMillan degree of H(s), which is 2.
    # D(s) is constructed to clear denominators and N(s) = D(s)H(s).
    # The key is to find a D(s) that is "minimal".
    
    # A minimal D(s) and the corresponding N(s) are:
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])
    
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("D(s) =")
    sympy.pretty_print(D)
    print("\nN(s) =")
    sympy.pretty_print(N)
    print("\n" + "="*40 + "\n")
    
    # Step 2: Verification
    print("Verification Step 1: Check if D(s)^-1 * N(s) == H(s)")
    # Calculate D_inv_N = D(s)^-1 * N(s)
    D_inv = D.inv()
    D_inv_N = D_inv * N
    
    # Simplify the resulting matrix
    D_inv_N_simplified = sympy.simplify(D_inv_N)
    
    print("D(s)^-1 * N(s) simplifies to:")
    sympy.pretty_print(D_inv_N_simplified)
    
    if D_inv_N_simplified == H:
        print("\nVerification successful: D(s)^-1 * N(s) equals H(s).")
    else:
        print("\nVerification failed: D(s)^-1 * N(s) does not equal H(s).")

    print("\n" + "="*40 + "\n")

    # Step 3: Coprimeness Check
    print("Verification Step 2: Check for left coprimeness.")
    print("This is done by checking if the GCD of all 2x2 minors of [D(s) N(s)] is a constant.")
    
    # Form the concatenated matrix [D(s) N(s)]
    DN_matrix = D.row_join(N)
    print("\n[D(s) N(s)] =")
    sympy.pretty_print(DN_matrix)
    
    # Calculate all 2x2 minors
    minors = []
    for j1 in range(DN_matrix.cols):
        for j2 in range(j1 + 1, DN_matrix.cols):
            sub_matrix = DN_matrix.col_slice(j1, j1+1).row_join(DN_matrix.col_slice(j2, j2+1))
            minors.append(sub_matrix.det())
            
    # Calculate the GCD of all minors
    if len(minors) > 1:
        gcd_of_minors = sympy.gcd(minors[0], minors[1])
        for i in range(2, len(minors)):
            gcd_of_minors = sympy.gcd(gcd_of_minors, minors[i])
    elif minors:
        gcd_of_minors = minors[0]
    else:
        gcd_of_minors = 0
        
    print(f"\nThe GCD of the 2x2 minors is: {gcd_of_minors}")

    # A non-zero constant GCD implies coprimeness
    if not gcd_of_minors.has(s) and gcd_of_minors != 0:
        print("Since the GCD is a non-zero constant, D(s) and N(s) are left coprime.")
        # Final answer format requires printing the equation elements
        print("\nFinal equation elements:")
        print("D(s) = ")
        for row in D.tolist():
            print(row)
        print("N(s) = ")
        for row in N.tolist():
            print(row)
    else:
        print("D(s) and N(s) are NOT left coprime.")

if __name__ == '__main__':
    solve_factorization()