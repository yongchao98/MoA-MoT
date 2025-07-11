import sympy as sp

def solve_coprime_factorization():
    """
    Calculates and prints a left coprime factorization of the transfer function H(s).
    The factorization is of the form H(s) = D(s)^-1 * N(s).
    """
    # Define the symbolic variable s
    s = sp.Symbol('s')

    # Define the transfer function matrix H(s)
    # H = sp.Matrix([
    #     [(s - 1) / (s + 1), 1],
    #     [2 / (s**2 - 1), 0]
    # ])

    # Step 1: Choose a constant non-singular matrix N(s).
    # This choice ensures the resulting factorization is coprime.
    N = sp.Matrix([
        [0, 2],
        [1, 1]
    ])

    # Step 2: From the relation D(s)H(s) = N(s), we derive expressions for the elements of D(s).
    # For the given H(s), we have n12 = d11 and n22 = d21.
    n11, n12 = N[0, 0], N[0, 1]
    n21, n22 = N[1, 0], N[1, 1]

    d11 = n12
    d21 = n22

    # Step 3: Solve for the remaining elements of D(s).
    # d12 = (n11(s^2-1) - d11(s-1)^2) / 2
    # d22 = (n21(s^2-1) - d21(s-1)^2) / 2
    d12_expr = (n11 * (s**2 - 1) - d11 * (s - 1)**2) / 2
    d22_expr = (n21 * (s**2 - 1) - d21 * (s - 1)**2) / 2
    
    d12 = sp.simplify(d12_expr)
    d22 = sp.simplify(d22_expr)

    # Step 4: Construct the D(s) matrix.
    D = sp.Matrix([
        [d11, d12],
        [d21, d22]
    ])

    # Step 5: Print the final factorization.
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("")
    
    # Pretty print D(s)
    d11_str = str(sp.expand(D[0,0]))
    d12_str = str(sp.expand(D[0,1]))
    d21_str = str(sp.expand(D[1,0]))
    d22_str = str(sp.expand(D[1,1]))
    
    print("D(s) =")
    print(f"[[ {d11_str:<18} , {d12_str:<18} ]]")
    print(f" [ {d21_str:<18} , {d22_str:<18} ]]")
    print("")

    # Pretty print N(s)
    n11_str = str(sp.expand(N[0,0]))
    n12_str = str(sp.expand(N[0,1]))
    n21_str = str(sp.expand(N[1,0]))
    n22_str = str(sp.expand(N[1,1]))

    print("N(s) =")
    print(f"[[ {n11_str:<18} , {n12_str:<18} ]]")
    print(f" [ {n21_str:<18} , {n22_str:<18} ]]")

solve_coprime_factorization()