import sympy

def coprime_factorization():
    """
    This function defines and prints the matrices for a left coprime factorization
    of the given transfer function H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # The resulting matrices D(s) and N(s) for the left coprime factorization
    # H(s) = D(s)^-1 * N(s)
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])

    # Print the final answer in a structured format
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")

    print("\nD(s) =")
    sympy.pprint(D)

    print("\nN(s) =")
    sympy.pprint(N)
    
    # Print each element of the final equation as requested
    print("\nIn element-wise form, the components of D(s) and N(s) are:")
    print(f"D(1,1) = {D[0,0]}")
    print(f"D(1,2) = {D[0,1]}")
    print(f"D(2,1) = {D[1,0]}")
    print(f"D(2,2) = {D[1,1]}")

    print(f"\nN(1,1) = {N[0,0]}")
    print(f"N(1,2) = {N[0,1]}")
    print(f"N(2,1) = {N[1,0]}")
    print(f"N(2,2) = {N[1,1]}")

coprime_factorization()
<<<D(s) = [[1, s - 1], [0, s**2 - 1]], N(s) = [[1, 1], [2, 0]]>>>