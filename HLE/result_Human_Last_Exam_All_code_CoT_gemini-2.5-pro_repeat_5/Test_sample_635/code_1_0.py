import sympy

def calculate_left_coprime_factorization():
    """
    This function defines and prints the matrices D(s) and N(s) for the
    left coprime factorization H(s) = D(s)^-1 * N(s).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # The derived left coprime factorization matrices are D(s) and N(s).
    # D(s) = [[1, s-1], [0, s^2-1]]
    # N(s) = [[1, 1], [2, 0]]
    D_s = sympy.Matrix([[1, s - 1],
                      [0, s**2 - 1]])

    N_s = sympy.Matrix([[1, 1],
                      [2, 0]])

    print("A left coprime factorization of H(s) is H(s) = D(s)^-1 * N(s), where:")
    print("\n")

    # The use_unicode=True flag ensures a nice print format in supporting terminals.
    print("D(s) =")
    sympy.pprint(D_s, use_unicode=True)
    print("\n")

    print("N(s) =")
    sympy.pprint(N_s, use_unicode=True)
    print("\n")
    
    # Print each element of the matrices as requested.
    print("--- The elements of the matrices are ---")
    print(f"D(1,1) = {D_s[0,0]}")
    print(f"D(1,2) = {D_s[0,1]}")
    print(f"D(2,1) = {D_s[1,0]}")
    print(f"D(2,2) = {D_s[1,1]}")
    print("")
    print(f"N(1,1) = {N_s[0,0]}")
    print(f"N(1,2) = {N_s[0,1]}")
    print(f"N(2,1) = {N_s[1,0]}")
    print(f"N(2,2) = {N_s[1,1]}")

calculate_left_coprime_factorization()