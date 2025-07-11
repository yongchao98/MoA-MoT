import sympy

def solve_factorization():
    """
    Calculates and prints a left coprime factorization for the given transfer function.
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # Define the calculated polynomial matrices D(s) and N(s)
    D = sympy.Matrix([[1, s - 1],
                      [0, s**2 - 1]])

    N = sympy.Matrix([[1, 1],
                      [2, 0]])

    # Print the factorization equation and the matrices
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    
    print("\nD(s) =")
    sympy.pprint(D, use_unicode=True)
    
    print("\nN(s) =")
    sympy.pprint(N, use_unicode=True)

    # The final equation is H(s) = D(s)^-1 * N(s), where each matrix is:
    print("\nIn equation form:")
    print("H(s) = ")
    sympy.pprint(D.inv(), use_unicode=True)
    print(" * ")
    sympy.pprint(N, use_unicode=True)


solve_factorization()