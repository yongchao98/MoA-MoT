import sympy

def solve_factorization():
    """
    This function calculates and displays a left coprime factorization for the given
    transfer function H(s).

    H(s) = [[(s-1)/(s+1), 1], [2/(s**2-1), 0]]
    We want to find D(s) and N(s) such that H(s) = D(s)^-1 * N(s)
    and D(s), N(s) are left coprime.

    Step 1: Find a common denominator factorization.
    The least common denominator is s^2 - 1.
    H(s) = (1/(s^2-1)) * [[(s-1)^2, s^2-1], [2, 0]]
    This gives an initial (non-coprime) factorization:
    D1(s) = [[s^2-1, 0], [0, s^2-1]]
    N1(s) = [[s^2-2*s+1, s^2-1], [2, 0]]

    Step 2: Find and cancel the Greatest Common Left Divisor (GCLD).
    Through analysis (e.g., using row reduction on [D1(s) N1(s)]), the GCLD G(s)
    can be found. The final coprime matrices are D(s) and N(s) where
    D1(s) = G(s)D(s) and N1(s) = G(s)N(s).

    The calculated coprime matrices are:
    D(s) = [[1, s-1], [0, s^2-1]]
    N(s) = [[1, 1], [2, 0]]
    """

    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # The resulting left coprime polynomial matrices D(s) and N(s)
    D_s = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N_s = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])

    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    
    print("\nD(s) =")
    sympy.pprint(D_s)
    
    print("\nN(s) =")
    sympy.pprint(N_s)

    # Optional: Verification step
    # H_calc = sympy.simplify(D_s.inv() * N_s)
    # print("\nVerification: D(s)^-1 * N(s) simplifies to:")
    # sympy.pprint(H_calc)
    # H_original = sympy.Matrix([[(s-1)/(s+1), 1], [2/(s**2-1), 0]])
    # print("\nWhich matches the original H(s):")
    # sympy.pprint(H_original)


if __name__ == "__main__":
    solve_factorization()
