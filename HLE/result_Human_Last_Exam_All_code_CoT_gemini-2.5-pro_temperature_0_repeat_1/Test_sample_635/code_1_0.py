import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and displays a left coprime factorization for the given transfer function.
    H(s) = D(s)^-1 * N(s)
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # The determined left coprime factorization matrices D(s) and N(s)
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])

    # Print the matrices D(s) and N(s) that form the factorization
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("\nD(s) =")
    sympy.pprint(D)
    
    # Output each element of D(s) for the equation
    print("\nEquation for D(s):")
    print(f"D(1,1) = {D[0,0]}")
    print(f"D(1,2) = {D[0,1]}")
    print(f"D(2,1) = {D[1,0]}")
    print(f"D(2,2) = {D[1,1]}")

    print("\nN(s) =")
    sympy.pprint(N)

    # Output each element of N(s) for the equation
    print("\nEquation for N(s):")
    print(f"N(1,1) = {N[0,0]}")
    print(f"N(1,2) = {N[0,1]}")
    print(f"N(2,1) = {N[1,0]}")
    print(f"N(2,2) = {N[1,1]}")

    # Verification (optional, but good practice)
    # H = D^-1 * N
    H_reconstructed = sympy.simplify(D.inv() * N)
    
    print("\nVerification:")
    print("H(s) reconstructed from D(s)^-1 * N(s):")
    sympy.pprint(H_reconstructed)

calculate_left_coprime_factorization()