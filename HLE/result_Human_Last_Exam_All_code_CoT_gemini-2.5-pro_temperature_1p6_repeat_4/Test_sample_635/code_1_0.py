import sympy

def solve_coprime_factorization():
    """
    Calculates and prints a left coprime factorization for the given transfer function.
    H(s) = D(s)^-1 * N(s)
    """
    s = sympy.symbols('s')

    # The transfer function matrix H(s)
    # H = sympy.Matrix([
    #     [(s - 1) / (s + 1), 1],
    #     [2 / (s**2 - 1), 0]
    # ])

    # Based on the derivation, the left coprime factorization H(s) = D(s)^-1 * N(s) is found.
    # We define the polynomial matrices D(s) and N(s).
    D = sympy.Matrix([
        [s + 1, 0],
        [1, s - 1]
    ])

    N = sympy.Matrix([
        [s - 1, s + 1],
        [1, 1]
    ])

    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("\nD(s) =")
    sympy.pprint(D)
    
    # We can represent the matrix as a list of lists for clear output of each element
    # as requested.
    d_list = D.tolist()
    print("\nEach number/expression in D(s):")
    for row in d_list:
        print(row)

    print("\nN(s) =")
    sympy.pprint(N)
    
    n_list = N.tolist()
    print("\nEach number/expression in N(s):")
    for row in n_list:
        print(row)

    # Verification (optional, for confidence)
    # H_calc = D.inv() * N
    # H_calc_simplified = sympy.simplify(H_calc)
    # print("\nVerification: D(s)^-1 * N(s) = ")
    # sympy.pprint(H_calc_simplified)

if __name__ == '__main__':
    solve_coprime_factorization()