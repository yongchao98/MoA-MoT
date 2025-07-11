import sympy

def solve_coprime_factorization():
    """
    Calculates a left coprime factorization H(s) = D(s)^-1 * N(s) for the given transfer function matrix.
    The method involves finding a suitable polynomial matrix D(s) that makes N(s) = D(s)H(s) a polynomial matrix,
    and then verifying the coprimeness of the pair (D(s), N(s)).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    print("The given transfer function is H(s):")
    sympy.pprint(H)
    print("\nWe are looking for a factorization H(s) = D(s)^-1 * N(s).\n")

    # From the condition D(s)H(s) = N(s), where D and N are polynomial matrices,
    # we derive conditions on the elements of D(s).
    # A systematic derivation leads to a simple choice for D(s) that satisfies the requirements.
    # The simplest non-trivial solution is found to be:
    
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    # Now, we calculate the corresponding N(s) matrix using N(s) = D(s)H(s)
    # sympy.simplify() is used to cancel out terms and show N(s) is polynomial.
    N = sympy.simplify(D * H)

    # The resulting matrices for the left coprime factorization are:
    print("A valid left coprime factorization is given by:")
    print("D(s) =")
    sympy.pprint(D)
    print("\nN(s) =")
    sympy.pprint(N)
    print("\n")

    # The final equation is:
    print("H(s) = D(s)^-1 * N(s) = ")
    # Print D inverse
    print("⎡ 1   s - 1⎤⁻¹   ⎡1  1⎤")
    print("⎢          ⎥      ⎢    ⎥")
    print("⎣ 0  s² - 1⎦      ⎣2  0⎦")


    # Verification steps
    # 1. Check if D(s)^-1 * N(s) equals H(s)
    H_reconstructed = sympy.simplify(D.inv() * N)
    is_correct_factorization = (H_reconstructed == H)
    
    # 2. Check if D(s) and N(s) are left coprime
    M = D.row_join(N)
    # The rank of M is 2 for all s if at least one 2x2 minor is a non-zero constant.
    minor_13 = M[:, [0, 2]].det() # Minor from columns 1 and 3 of [D N]
    is_coprime = minor_13 != 0 and minor_13.is_constant()

    print("\n--- Verification ---")
    print(f"Factorization Correct: D(s)^-1 * N(s) = H(s)?  -> {is_correct_factorization}")
    print(f"D(s) and N(s) are Left Coprime? -> {is_coprime}")
    if is_coprime:
        print("The rank of [D(s) N(s)] is 2 for all s. For example, the minor formed by columns 1 and 3 has a determinant of 2.")

solve_coprime_factorization()