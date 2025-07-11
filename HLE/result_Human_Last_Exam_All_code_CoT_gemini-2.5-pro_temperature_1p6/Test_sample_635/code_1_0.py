import sympy

def solve_left_coprime_factorization():
    """
    Calculates and verifies a left coprime factorization for the given
    transfer function H(s).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # The transfer function H(s) given in the problem
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # Our calculated left coprime factorization matrices D(s) and N(s)
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])

    # --- Output the Results ---
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    
    # Pretty print the final D(s) matrix
    print("\nD(s) =")
    sympy.pretty_print(D)
    print(f"\nD(s)[0,0] = {D[0,0]}")
    print(f"D(s)[0,1] = {D[0,1]}")
    print(f"D(s)[1,0] = {D[1,0]}")
    print(f"D(s)[1,1] = {D[1,1]}")

    # Pretty print the final N(s) matrix
    print("\nN(s) =")
    sympy.pretty_print(N)
    print(f"\nN(s)[0,0] = {N[0,0]}")
    print(f"N(s)[0,1] = {N[0,1]}")
    print(f"N(s)[1,0] = {N[1,0]}")
    print(f"N(s)[1,1] = {N[1,1]}")

    # --- Verification (Optional) ---
    print("\n\n--- Verification Steps ---")

    # 1. Verify D(s) * H(s) = N(s)
    print("\n1. Verifying D(s) * H(s) = N(s):")
    calculated_N = sympy.simplify(D * H)
    is_factorization_correct = (calculated_N == N)
    print("Result of D(s) * H(s):")
    sympy.pretty_print(calculated_N)
    print(f"Is the factorization D(s)H(s) = N(s) correct? {is_factorization_correct}")

    # 2. Verify D(s) and N(s) are left coprime
    print("\n2. Verifying D(s) and N(s) are left coprime:")
    M = sympy.BlockMatrix([[D, N]]).as_explicit()
    print("The matrix [D(s) N(s)] is:")
    sympy.pretty_print(M)
    
    # Calculate all 2x2 minors of M
    minors = []
    for c1 in range(M.cols):
        for c2 in range(c1 + 1, M.cols):
            minor_det = sympy.simplify(M[:, [c1, c2]].det())
            minors.append(minor_det)
    
    # The GCD of the minors must be a non-zero constant for coprimeness
    gcd_of_minors = sympy.gcd_list(minors)
    is_coprime = gcd_of_minors.is_constant and not gcd_of_minors.is_zero
    
    print(f"\nThe greatest common divisor (GCD) of all 2x2 minors is: {gcd_of_minors}")
    if is_coprime:
        print("Since the GCD is a non-zero constant, D(s) and N(s) are left coprime.")
    else:
        print("Since the GCD is not a non-zero constant, D(s) and N(s) are not left coprime.")

if __name__ == '__main__':
    solve_left_coprime_factorization()
