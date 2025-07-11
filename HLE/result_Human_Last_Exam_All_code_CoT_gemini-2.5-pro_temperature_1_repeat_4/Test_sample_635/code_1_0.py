import sympy

def calculate_left_coprime_factorization():
    """
    Calculates a left coprime factorization of the transfer function H(s).
    The factorization is of the form H(s) = D(s)^-1 * N(s).
    """
    # Initialize the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    print("Given the transfer function H(s):")
    sympy.pprint(H, use_unicode=False)
    print("\nWe want to find a left coprime factorization H(s) = D(s)^-1 * N(s).\n")

    # Step 1: Analyze H(s) for unstable poles.
    # The denominators are (s+1) and (s^2 - 1) = (s-1)(s+1).
    # The poles are at s = -1 (stable) and s = 1 (unstable).
    # The unstable pole at s=1 appears in the H[1,0] element.

    # Step 2: Construct D(s) to cancel the instability.
    # The first row of H(s) is composed of stable and proper functions, so we can
    # choose the first row of D(s) to be [1, 0] which leaves the first row of N(s) stable.
    # The second row of H(s) contains the unstable pole. We need to design the second
    # row of D(s) to cancel it. We choose a diagonal D(s) for simplicity.
    # D(s) = diag(d11, d22)
    # d11(s) = 1.
    # N[1,0] = d22(s) * H[1,0] = d22(s) * 2/((s-1)(s+1)).
    # To cancel the unstable pole at s=1, d22(s) must have a zero at s=1.
    # A simple, stable, and proper choice for d22(s) is (s-1)/(s+a) with a>0.
    # We choose a=1 for simplicity.
    
    d22 = (s - 1) / (s + 1)

    # Form the D(s) matrix
    D = sympy.Matrix([
        [1, 0],
        [0, d22]
    ])

    # Step 3: Calculate N(s) = D(s) * H(s)
    N = D * H

    # Simplify the resulting matrix N(s)
    N = sympy.simplify(N)

    # Step 4: Print the final results in the required format.
    print("A left coprime factorization is found with D(s) and N(s) as follows:")
    print("------------------------------------------------------------------")
    
    print("D(s) =")
    sympy.pprint(D, use_unicode=False)
    print("\nN(s) =")
    sympy.pprint(N, use_unicode=False)
    print("\n------------------------------------------------------------------")
    print("The final equation is H(s) = D(s)^-1 * N(s), where the elements of each matrix are:\n")

    print(f"D(s)_11 = {D[0,0]}")
    print(f"D(s)_12 = {D[0,1]}")
    print(f"D(s)_21 = {D[1,0]}")
    print(f"D(s)_22 = {D[1,1]}")
    print("")
    print(f"N(s)_11 = {N[0,0]}")
    print(f"N(s)_12 = {N[0,1]}")
    print(f"N(s)_21 = {N[1,0]}")
    print(f"N(s)_22 = {N[1,1]}")


if __name__ == '__main__':
    calculate_left_coprime_factorization()
