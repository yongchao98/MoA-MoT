import sympy

def calculate_left_coprime_factorization():
    """
    This function calculates a left coprime factorization H(s) = D(s)^-1 * N(s)
    for the given transfer function matrix H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    # H(s) = [[ (s-1)/(s+1),   1 ],
    #         [ 2/(s^2-1),     0 ]]
    H = sympy.Matrix([
        [(s-1)/(s+1), 1],
        [2/(s**2-1), 0]
    ])

    print("Calculating a left coprime factorization for H(s):")
    sympy.pprint(H)

    # To make D(s) * H(s) a polynomial matrix, we need to choose D(s) to cancel
    # the denominators in H(s).
    # Let D(s) = [[d11, d12], [d21, d22]].
    # N(s) = D(s) * H(s)
    # n11 = d11*(s-1)/(s+1) + d12*2/(s^2-1)
    # n21 = d21*(s-1)/(s+1) + d22*2/(s^2-1)
    # We choose a simple D(s) that satisfies the polynomial condition for N(s).
    # A simple choice that works is:
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    # Calculate N(s) = D(s) * H(s) and simplify
    N = sympy.simplify(D * H)

    print("\nA left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")

    # Print D(s) and N(s)
    print("\nD(s) =")
    sympy.pprint(D)

    print("\nN(s) =")
    sympy.pprint(N)

    # --- Verification ---
    print("\n--- Verification Steps ---")
    # 1. Verify that D(s)^-1 * N(s) = H(s)
    H_verified = sympy.simplify(D.inv() * N)
    print("1. Verifying the factorization D(s)^-1 * N(s):")
    sympy.pprint(H_verified)
    is_correct = (H_verified == H)
    print(f"The factorization is correct: {is_correct}")

    # 2. Verify left coprimeness of D(s) and N(s)
    # Form the augmented matrix M = [D(s) N(s)]
    M = D.row_join(N)
    print("\n2. Verifying left coprimeness by checking the rank of [D(s) N(s)]:")
    sympy.pprint(M)
    
    # Calculate the 2x2 minors of M. If any minor is a non-zero constant,
    # the matrix has full row rank for all s, and D, N are left coprime.
    minor_1_3 = M[:, [0, 2]].det() # Minor from column 1 and 3
    print(f"\nA 2x2 minor of this matrix is det([[D_11, N_11], [D_21, N_21]]) = {minor_1_3}.")

    if minor_1_3.is_constant() and minor_1_3 != 0:
        print("Since this minor is a non-zero constant, [D(s) N(s)] has full rank for all s.")
        print("Therefore, D(s) and N(s) are left coprime.")
    else:
        print("The matrices D(s) and N(s) are not left coprime.")
    
    # Print final equation with each number explicitly
    print("\n--- Final Equation ---")
    print("The factorization is given by H(s) = D(s)^-1 * N(s), where:")
    print("\nD(s) = ")
    print(f"[[ {D[0,0]} , {D[0,1]} ]]")
    print(f"[ [ {D[1,0]} , {D[1,1]} ]]")

    print("\nN(s) = ")
    print(f"[[ {N[0,0]} , {N[0,1]} ]]")
    print(f"[ [ {N[1,0]} , {N[1,1]} ]]")


if __name__ == "__main__":
    calculate_left_coprime_factorization()
