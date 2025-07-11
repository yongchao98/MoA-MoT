import sympy

def calculate_left_coprime_factorization():
    """
    Calculates a left coprime factorization H(s) = D(s)^-1 * N(s) for the given transfer function matrix.
    """
    # Define the symbolic variable s
    s = sympy.symbols('s')

    # Define the transfer function H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # Step 1: Find an initial, non-coprime factorization D1(s), N1(s)
    # This is done by finding a common denominator for each row.
    d1 = s + 1
    d2 = s**2 - 1
    D1 = sympy.Matrix([
        [d1, 0],
        [0, d2]
    ])
    N1 = sympy.simplify(D1 * H)

    # Step 2: Check for coprimeness by examining the rank at the roots of det(D1(s))
    # det(D1(s)) = (s+1)*(s^2-1), roots are s=1, s=-1.
    M1 = D1.row_join(N1)
    
    # At s = -1, the matrix M1 loses rank
    M1_at_minus_1 = M1.subs(s, -1)
    
    # Step 3: Extract the common left divisor
    # Find the left null space of M1 at s = -1
    # The nullspace of M.T gives the left nullspace of M
    left_nullspace_basis = M1_at_minus_1.T.nullspace()
    
    # The basis vector is [1, 1]. We form a unimodular transformation matrix U.
    v = left_nullspace_basis[0].T  # This gives Matrix([[1, 1]])
    U = sympy.eye(2)
    U[1, :] = v # Replace the second row with v

    # Apply the transformation to D1 and N1
    D2 = sympy.simplify(U * D1)
    N2 = sympy.simplify(U * N1)

    # The second row is now divisible by (s+1). We extract this factor.
    D_final = D2.copy()
    N_final = N2.copy()
    
    # Divide the second row by (s+1)
    row_to_divide = D2.row(1)
    D_final[1, :] = row_to_divide.applyfunc(lambda x: sympy.cancel(x / (s + 1)))
    
    row_to_divide_N = N2.row(1)
    N_final[1, :] = row_to_divide_N.applyfunc(lambda x: sympy.cancel(x / (s + 1)))
    
    # The final matrices are D_final and N_final
    D = D_final
    N = N_final

    print("The left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("\nD(s) =")
    sympy.pprint(D)
    print("\nN(s) =")
    sympy.pprint(N)
    
    # Verification of the result
    # H_check = D.inv() * N
    # assert sympy.simplify(H_check - H) == sympy.zeros(2,2)

if __name__ == '__main__':
    calculate_left_coprime_factorization()