import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and verifies a left coprime factorization of the given
    transfer function matrix H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.symbols('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    print("Original Transfer Function H(s):")
    sympy.pretty_print(H)
    print("-" * 30)

    # The left coprime factorization H(s) = D(s)^-1 * N(s) is found by
    # clearing denominators and removing common left divisors.
    # The resulting matrices are D(s) and N(s).

    # Define the calculated denominator matrix D(s)
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    # Define the calculated numerator matrix N(s)
    N = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])

    print("Calculated Left Coprime Factorization:")
    print("D(s) =")
    sympy.pretty_print(D)
    print("\nN(s) =")
    sympy.pretty_print(N)
    print("-" * 30)

    # --- Verification Step 1: Check if D^-1 * N = H ---
    print("Verification 1: Check if D(s)^-1 * N(s) == H(s)")
    # Calculate D^-1 * N
    H_reconstructed = sympy.simplify(D.inv() * N)
    if H_reconstructed == H:
        print("Result: Success! D(s)^-1 * N(s) is equal to H(s).")
    else:
        print("Result: Failure! D(s)^-1 * N(s) is NOT equal to H(s).")
        print("D(s)^-1 * N(s) evaluates to:")
        sympy.pretty_print(H_reconstructed)
    print("-" * 30)
    
    # --- Verification Step 2: Check for left coprimeness ---
    # The pair (D, N) is left coprime if the matrix [D(s) N(s)] has
    # full row rank for all s. We only need to check this at the roots
    # of det(D(s)).
    print("Verification 2: Check for coprimeness")
    det_D = D.det()
    print(f"det(D(s)) = {det_D}")

    # Find the roots of the determinant
    roots = sympy.solve(det_D, s)
    print(f"Roots of det(D(s)) are: {roots}")

    # Form the composite matrix [D(s) N(s)]
    M = D.row_join(N)
    print("\nComposite matrix [D(s) N(s)]:")
    sympy.pretty_print(M)

    is_coprime = True
    for root in roots:
        M_at_root = M.subs(s, root)
        rank_at_root = M_at_root.rank()
        print(f"\nAt s = {root}:")
        print("Matrix [D(s) N(s)] becomes:")
        sympy.pretty_print(M_at_root)
        print(f"Rank = {rank_at_root}")
        if rank_at_root < D.shape[0]:
            is_coprime = False
            print(f"Rank is less than {D.shape[0]}. D(s) and N(s) are NOT coprime.")
            
    if is_coprime:
        print("\nResult: Success! The rank is full at all roots. D(s) and N(s) are left coprime.")
    print("-" * 30)

    # --- Final Result ---
    print("The final factorization is H(s) = D(s)^-1 * N(s):")
    # Print the equation with each "number" (matrix element)
    print(f"\nH(s) = \n")
    sympy.pretty_print(D)
    print("     ^-1")
    print("      *")
    sympy.pretty_print(N)


if __name__ == '__main__':
    calculate_left_coprime_factorization()
<<<D(s) = Matrix([[1, s - 1], [0, s**2 - 1]]), N(s) = Matrix([[1, 1], [2, 0]])>>>