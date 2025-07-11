import sympy

def solve_left_coprime_factorization():
    """
    Calculates and displays a left coprime factorization for the given transfer function H(s).
    The factorization is of the form H(s) = D(s)^-1 * N(s).
    """
    # Initialize the symbolic variable and printing
    s = sympy.symbols('s')
    sympy.init_printing(use_unicode=True, wrap_line=False)

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])
    print("Given transfer function H(s):")
    sympy.pprint(H)
    print("-" * 40)

    # Step 1: Find an initial factorization by clearing denominators row-by-row.
    # The least common denominator for row 1 is (s+1).
    # The least common denominator for row 2 is (s^2-1).
    D_prelim = sympy.Matrix([
        [s + 1, 0],
        [0, s**2 - 1]
    ])
    
    # Calculate N_prelim = D_prelim * H
    N_prelim = sympy.simplify(D_prelim * H)

    print("Step 1: Initial (non-coprime) factorization H(s) = D_prelim(s)^-1 * N_prelim(s)")
    print("D_prelim(s) =")
    sympy.pprint(D_prelim)
    print("\nN_prelim(s) =")
    sympy.pprint(N_prelim)
    print("-" * 40)

    # Step 2: Check for coprimeness at the roots of det(D_prelim(s)).
    # det(D_prelim(s)) = (s+1)*(s^2-1). The roots are s=1 and s=-1.
    M = D_prelim.row_join(N_prelim)
    
    # Check at s = -1
    M_at_minus_1 = M.subs(s, -1)
    rank = M_at_minus_1.rank()
    print("Step 2: Check coprimeness. A rank deficiency is found at s = -1:")
    print("[D_prelim(-1) N_prelim(-1)] =")
    sympy.pprint(M_at_minus_1)
    print(f"The rank is {rank}, which is less than 2. The factorization is not coprime.")
    print("The dependency is Row 1 + Row 2 = 0 at s = -1.")
    print("-" * 40)

    # Step 3: Remove the common factor (s+1) using a unimodular transformation.
    # The dependency R1+R2=0 suggests the transformation matrix U.
    U = sympy.Matrix([
        [1, 1],
        [0, 1]
    ])
    
    D_prime = sympy.simplify(U * D_prelim)
    N_prime = sympy.simplify(U * N_prelim)
    
    print("Step 3: Apply unimodular transformation to extract the common factor.")
    print("The first row of [U*D_prelim, U*N_prelim] now has a common factor of (s+1).")

    # Extract the common factor G(s) = diag(s+1, 1) from the left.
    G_inv = sympy.Matrix([
        [1 / (s + 1), 0],
        [0, 1]
    ])
    
    D_final = sympy.simplify(G_inv * D_prime)
    N_final = sympy.simplify(G_inv * N_prime)
    
    print("After extracting the factor, we get the final coprime matrices.")
    print("-" * 40)

    # Final Result
    print("Final Left Coprime Factorization H(s) = D(s)^-1 * N(s)")
    print("\nD(s) =")
    sympy.pprint(D_final)
    print("\nN(s) =")
    sympy.pprint(N_final)

    # Print the final equation with each term
    print("\nThe final equation is:")
    
    # Convert matrices to string lists for formatting
    d_str = [[str(e) for e in row] for row in D_final.tolist()]
    n_str = [[str(e) for e in row] for row in N_final.tolist()]
    
    # Calculate column widths for alignment
    d_widths = [max(len(d_str[i][j]) for i in range(2)) for j in range(2)]
    n_widths = [max(len(n_str[i][j]) for i in range(2)) for j in range(2)]

    # Print formatted equation
    print(f"⎡ {d_str[0][0]:<{d_widths[0]}}  {d_str[0][1]:<{d_widths[1]}} ⎤-1   ⎡ {n_str[0][0]:<{n_widths[0]}}  {n_str[0][1]:<{n_widths[1]}} ⎤")
    print(f"⎢ {' ' * d_widths[0]}  {' ' * d_widths[1]} ⎥   * ⎢ {' ' * n_widths[0]}  {' ' * n_widths[1]} ⎥")
    print(f"⎣ {d_str[1][0]:<{d_widths[0]}}  {d_str[1][1]:<{d_widths[1]}} ⎦      ⎣ {n_str[1][0]:<{n_widths[0]}}  {n_str[1][1]:<{n_widths[1]}} ⎦")

solve_left_coprime_factorization()