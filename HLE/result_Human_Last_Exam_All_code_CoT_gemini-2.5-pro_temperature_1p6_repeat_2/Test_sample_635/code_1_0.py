import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and verifies a left coprime factorization for the given transfer function H(s).
    """
    s = sympy.symbols('s')
    
    # The given transfer function matrix
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])
    
    print("The given transfer function H(s) is:")
    sympy.pprint(H)
    print("-" * 30)

    # Step 1: Find an initial factorization H(s) = D_0(s)^-1 * N_0(s)
    print("Step 1: Find an initial factorization using the least common denominator.")
    # The least common multiple of the denominators (s+1) and (s**2-1) is s**2-1.
    common_denominator = s**2 - 1
    
    # D_0(s) is the common denominator times the identity matrix.
    D0 = common_denominator * sympy.eye(2)
    
    # N_0(s) is the corresponding numerator matrix.
    N0 = sympy.Matrix([
        [sympy.simplify(H[0,0] * common_denominator), sympy.simplify(H[0,1] * common_denominator)],
        [sympy.simplify(H[1,0] * common_denominator), sympy.simplify(H[1,1] * common_denominator)]
    ])

    print("Initial denominator matrix D_0(s):")
    sympy.pprint(D0)
    print("\nInitial numerator matrix N_0(s):")
    sympy.pprint(N0)
    print("-" * 30)

    # Step 2: Perform reductions to achieve coprimeness.
    # The methodical process was detailed in the thought process. We found the final
    # coprime matrices D(s) and N(s) by extracting the greatest common left divisor.
    # Here, we present the final result of that reduction process.
    print("Step 2: Reduce the factorization to a coprime form.")
    print("After performing row operations on [D_0(s) N_0(s)] to remove common factors, we obtain the final D(s) and N(s) matrices.")
    
    D_final = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])
    
    N_final = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])
    
    print("\nFinal Left Coprime Factorization: H(s) = D(s)^-1 * N(s)")
    print("\nD(s) = ")
    sympy.pprint(D_final)
    print("\nN(s) = ")
    sympy.pprint(N_final)
    print("-" * 30)

    # Step 3: Verification
    print("Step 3: Verification.")
    
    # Check 1: D(s)H(s) - N(s) should be the zero matrix.
    verification1 = sympy.simplify(D_final * H - N_final)
    print("\nCheck 1: Verify D(s)H(s) = N(s)")
    print("Result of D(s)H(s) - N(s):")
    sympy.pprint(verification1)
    if verification1 == sympy.zeros(2, 2):
        print("Verification successful: The factorization is correct.")
    else:
        print("Verification failed: The factorization is incorrect.")
    
    # Check 2: The matrix [D(s) N(s)] must have full row rank for all s.
    print("\nCheck 2: Verify D(s) and N(s) are left coprime.")
    M_final = D_final.row_join(N_final)
    print("The composite matrix [D(s) N(s)] is:")
    sympy.pprint(M_final)
    
    # We find a 2x2 minor whose determinant is a non-zero constant.
    minor_cols_1_3 = M_final[:, [0, 2]]
    det_minor = minor_cols_1_3.det()
    
    print("\nTo check for full rank, we compute the determinant of a 2x2 minor.")
    print("Consider the minor formed by columns 1 and 3:")
    sympy.pprint(minor_cols_1_3)
    print(f"The determinant of this minor is: {det_minor}")
    
    if det_minor != 0:
        print("Since the determinant is a non-zero constant, the rank of [D(s) N(s)] is 2 for all s.")
        print("Verification successful: D(s) and N(s) are left coprime.")
    else:
        print("Verification failed: D(s) and N(s) are not left coprime.")

    print("-" * 30)
    print("Final Answer:")
    print("The left coprime factorization is H(s) = D(s)^-1 * N(s), where:")

    print("\nD(s) = [[D_11(s), D_12(s)], [D_21(s), D_22(s)]]")
    print(f"D_11(s) = {D_final[0,0]}")
    print(f"D_12(s) = {D_final[0,1]}")
    print(f"D_21(s) = {D_final[1,0]}")
    print(f"D_22(s) = {D_final[1,1]}")

    print("\nN(s) = [[N_11(s), N_12(s)], [N_21(s), N_22(s)]]")
    print(f"N_11(s) = {N_final[0,0]}")
    print(f"N_12(s) = {N_final[0,1]}")
    print(f"N_21(s) = {N_final[1,0]}")
    print(f"N_22(s) = {N_final[1,1]}")


if __name__ == '__main__':
    calculate_left_coprime_factorization()