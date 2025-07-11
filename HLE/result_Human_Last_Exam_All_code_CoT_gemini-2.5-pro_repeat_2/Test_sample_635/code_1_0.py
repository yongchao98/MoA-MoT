import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and prints a left coprime factorization H(s) = D(s)^-1 * N(s)
    for the given transfer function H(s).
    """
    # Define the symbol s
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # Step 1: Find an initial (not necessarily coprime) left factorization.
    # We choose D_0(s) as a diagonal matrix of row denominators.
    d1 = s + 1
    d2 = s**2 - 1
    D_0 = sympy.Matrix([[d1, 0], [0, d2]])

    # Calculate N_0(s) = D_0(s) * H(s)
    N_0 = sympy.simplify(D_0 * H)

    # Step 2: The initial factorization [D_0(s), N_0(s)] is not coprime.
    # A common factor exists related to (s+1). We use row operations to extract it.
    # At s = -1, row_2 of [D_0 N_0] is -1 times row_1.
    # The operation R2 -> R2 + R1 will reveal the common factor.
    # This is done by left-multiplying with a unimodular matrix E.
    E = sympy.Matrix([[1, 0], [1, 1]])
    D_1 = sympy.simplify(E * D_0)
    N_1 = sympy.simplify(E * N_0)

    # The second row of D_1 and N_1 is now divisible by (s+1).
    # We extract this common left divisor G(s) = diag(1, s+1).
    # The final coprime matrices are D = G^-1 * D_1 and N = G^-1 * N_1.
    G_inv = sympy.Matrix([[1, 0], [0, 1/(s+1)]])
    
    D = sympy.simplify(G_inv * D_1)
    N = sympy.simplify(G_inv * N_1)

    # Step 3: Print the resulting coprime matrices D(s) and N(s).
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")

    print("\nD(s) = ")
    sympy.pprint(D)

    print("\nN(s) = ")
    sympy.pprint(N)

    # Output each element of the matrices explicitly as requested.
    print("\nThe elements of the matrices are:")
    print(f"D_11(s) = {D[0,0]}")
    print(f"D_12(s) = {D[0,1]}")
    print(f"D_21(s) = {D[1,0]}")
    print(f"D_22(s) = {D[1,1]}")
    print("")
    print(f"N_11(s) = {N[0,0]}")
    print(f"N_12(s) = {N[0,1]}")
    print(f"N_21(s) = {N[1,0]}")
    print(f"N_22(s) = {N[1,1]}")

calculate_left_coprime_factorization()