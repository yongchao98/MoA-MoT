import sympy

def solve_factorization():
    """
    Calculates and prints a left coprime factorization for the given transfer function.
    """
    s = sympy.Symbol('s')

    # The factorization is found to be H(s) = D(s)^-1 * N(s)
    # where D(s) and N(s) are the polynomial matrices defined below.

    D_s = sympy.Matrix([
        [1, s - 1],
        [s + 1, 0]
    ])

    N_s = sympy.Matrix([
        [1, 1],
        [s - 1, s + 1]
    ])
    
    # We present the final factorization H(s) = D(s)^-1 * N(s) by printing D(s) and N(s)
    # in a clear, formatted way.
    
    print("A left coprime factorization is H(s) = D(s)^-1 * N(s), where:")
    
    # The requirement is to output each number/expression in the final equation.
    # The "final equation" defines D(s) and N(s).
    
    print("\nD(s) =")
    print(f"[{D_s[0, 0]:>7} {D_s[0, 1]:>7}]")
    print(f"[{D_s[1, 0]:>7} {D_s[1, 1]:>7}]")
    
    print("\nN(s) =")
    print(f"[{N_s[0, 0]:>7} {N_s[0, 1]:>7}]")
    print(f"[{N_s[1, 0]:>7} {N_s[1, 1]:>7}]")


solve_factorization()