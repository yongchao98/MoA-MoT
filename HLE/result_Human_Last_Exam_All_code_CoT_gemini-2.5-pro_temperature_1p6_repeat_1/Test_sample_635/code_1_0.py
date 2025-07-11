import sympy

def solve_factorization():
    """
    Calculates and prints a left coprime factorization of the given transfer function.
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the matrices D(s) and N(s) based on the derivation.
    # We found two independent row solutions for D(s): [s+1, 0] and [1, s-1].
    # Let's use these to construct D(s).
    D_s = sympy.Matrix([
        [s + 1, 0],
        [1, s - 1]
    ])

    # Corresponding N(s) is calculated as:
    # Row 1 of N: n11 = s-1, n12 = s+1
    # Row 2 of N: n21 = 1,   n22 = 1
    N_s = sympy.Matrix([
        [s - 1, s + 1],
        [1, 1]
    ])

    # Print the result in the specified format H(s) = D(s)^-1 * N(s)
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    print("")

    # Print D(s)
    print("D(s) = ")
    for i in range(D_s.rows):
        row_str = "[ "
        for j in range(D_s.cols):
            row_str += str(D_s[i, j])
            if j < D_s.cols - 1:
                row_str += ", "
        row_str += " ]"
        print(row_str)
    
    print("")

    # Print N(s)
    print("N(s) = ")
    for i in range(N_s.rows):
        row_str = "[ "
        for j in range(N_s.cols):
            row_str += str(N_s[i, j])
            if j < N_s.cols - 1:
                row_str += ", "
        row_str += " ]"
        print(row_str)

solve_factorization()