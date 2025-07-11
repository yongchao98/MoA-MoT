import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and prints a left coprime factorization H(s) = D(s)^-1 * N(s)
    for the given transfer function matrix H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # As derived in the explanation, a polynomial basis for the rows of D(s)
    # that clears the denominators in H(s) is {[s+1, 0], [1, s-1]}.
    # We construct D(s) using these basis vectors as its rows.
    row1 = sympy.Matrix([[s + 1, 0]])
    row2 = sympy.Matrix([[1, s - 1]])
    D = sympy.Matrix.vstack(row1, row2)

    # Calculate N(s) = D(s) * H(s).
    # The simplify() function is used to resolve the expressions and show
    # that the result is indeed a polynomial matrix.
    N = sympy.simplify(D * H)

    # --- Output the results ---
    print("The given transfer function H(s) is:")
    sympy.pprint(H, use_unicode=True)
    print("\nA left coprime factorization H(s) = D(s)^-1 * N(s) is found.")
    print("The matrices D(s) and N(s) are:\n")

    print("D(s) =")
    sympy.pprint(D, use_unicode=True)
    print("\nN(s) =")
    sympy.pprint(N, use_unicode=True)

    print("\nSo, the final equation is H(s) = D(s)^-1 * N(s) where:")
    
    # Printing each matrix element as requested.
    d11, d12 = D.row(0)
    d21, d22 = D.row(1)
    n11, n12 = N.row(0)
    n21, n22 = N.row(1)

    print(f"\nD(s) = [[{d11}, {d12}], [{d21}, {d22}]]")
    print(f"N(s) = [[{n11}, {n12}], [{n21}, {n22}]]")


if __name__ == '__main__':
    calculate_left_coprime_factorization()