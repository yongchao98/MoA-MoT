import sympy as sp

def solve_left_coprime_factorization():
    """
    Calculates and prints a left coprime factorization H(s) = D(s)^-1 * N(s)
    for the given transfer function.
    """
    # Define the symbolic variable s
    s = sp.symbols('s')

    # The given transfer function matrix H(s)
    H = sp.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # A left coprime factorization is found to be D(s) and N(s)
    # D(s) is constructed to cancel the denominators in H(s)
    D = sp.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    # N(s) is calculated as D(s) * H(s)
    N = sp.simplify(D * H)

    # --- Verification (optional, for confidence) ---
    # 1. Check if factorization is correct: H = D^-1 * N
    is_correct = (sp.simplify(D.inv() * N) == H)
    # 2. Check for coprimeness: GCD of all 2x2 minors of [D N] must be a constant.
    M = D.row_join(N)
    minors = M.minors(M.shape[0])
    gcd_of_minors = sp.gcd_list(minors)
    is_coprime = gcd_of_minors.is_constant() and gcd_of_minors != 0

    print("The factorization H(s) = D(s)^-1 * N(s) is:")
    if is_correct and is_coprime:
        print("Status: Verified as correct and left coprime.")
    else:
        print("Status: Verification failed.")
    print("-" * 40)
    
    # --- Final Answer ---
    print("The denominator polynomial matrix D(s) is:")
    sp.pprint(D)
    print("\nThe numerator polynomial matrix N(s) is:")
    sp.pprint(N)
    
    print("\n--- Final Equation Components ---")
    print("\nD(s) matrix elements:")
    print(f"D_11(s) = {D[0,0]}")
    print(f"D_12(s) = {D[0,1]}")
    print(f"D_21(s) = {D[1,0]}")
    print(f"D_22(s) = {D[1,1]}")

    print("\nN(s) matrix elements:")
    print(f"N_11(s) = {N[0,0]}")
    print(f"N_12(s) = {N[0,1]}")
    print(f"N_21(s) = {N[1,0]}")
    print(f"N_22(s) = {N[1,1]}")

if __name__ == '__main__':
    solve_left_coprime_factorization()
