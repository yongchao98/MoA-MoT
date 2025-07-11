import sympy

def calculate_left_coprime_factorization():
    """
    This function calculates a left coprime factorization H(s) = D(s)^-1 * N(s)
    for the given transfer function matrix H(s).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s-1)/(s+1), 1],
        [2/(s**2-1), 0]
    ])

    print("Step 1: Find an initial, non-coprime factorization D0(s) and N0(s).")
    # Choose D0(s) as a diagonal matrix of row-wise LCMs of denominators
    d1 = sympy.lcm(sympy.denom((s-1)/(s+1)), sympy.denom(1))
    d2 = sympy.lcm(sympy.denom(2/(s**2-1)), sympy.denom(0))
    D0 = sympy.diag(d1, d2)
    
    # Calculate N0(s) = D0(s) * H(s)
    N0 = sympy.simplify(D0 * H)

    print("Initial D0(s):")
    sympy.pprint(D0)
    print("\nInitial N0(s):")
    sympy.pprint(N0)
    print("-" * 30)

    print("Step 2: Check for coprimeness by testing the rank of [D0(s) | N0(s)].")
    M0 = D0.row_join(N0)
    det_D0 = D0.det()
    problematic_s_roots = sympy.roots(det_D0, s)
    
    print(f"The rank can drop at the roots of det(D0(s)): {list(problematic_s_roots.keys())}")
    
    # Check the rank at the problematic root s = -1
    s_val = -1
    M0_at_s_val = M0.subs(s, s_val)
    rank_at_s_val = M0_at_s_val.rank()
    
    print(f"At s = {s_val}, the rank of [D0(s) | N0(s)] is {rank_at_s_val}.")
    if rank_at_s_val < M0.rows:
        print("The rank is less than the number of rows, so the factorization is not coprime.")
    print("-" * 30)

    print("Step 3: Extract the common left divisor to find the coprime factorization.")
    # At s = -1, we see that R1 + R2 = 0. This implies R1(s) + R2(s) has a factor of (s+1).
    # We form a new second row by adding the first and second rows.
    R1 = M0.row(0)
    R2 = M0.row(1)
    new_R2 = R1 + R2

    # We divide the new second row by the common factor (s+1)
    factored_R2 = sympy.Matrix([sympy.simplify(elem / (s + 1)) for elem in new_R2]).T

    # The new coprime matrices D(s) and N(s) are formed from the old R1 and the new factored R2
    D = R1[:, 0:2].col_join(factored_R2[:, 0:2])
    N = R1[:, 2:4].col_join(factored_R2[:, 2:4])

    print("Final Result:")
    print("D(s) =")
    sympy.pprint(D)
    print("\nN(s) =")
    sympy.pprint(N)
    print("-" * 30)

    print("Final Equation H(s) = D(s)^-1 * N(s):")
    
    # Generate pretty string representations of the matrices
    h_str = f"[[{H[0,0]}, {H[0,1]}], [{H[1,0]}, {H[1,1]}]]"
    d_str = f"[[{D[0,0]}, {D[0,1]}], [{D[1,0]}, {D[1,1]}]]"
    n_str = f"[[{N[0,0]}, {N[0,1]}], [{N[1,0]}, {N[1,1]}]]"

    print(f"{h_str} = {d_str}^-1 * {n_str}")


if __name__ == '__main__':
    calculate_left_coprime_factorization()
