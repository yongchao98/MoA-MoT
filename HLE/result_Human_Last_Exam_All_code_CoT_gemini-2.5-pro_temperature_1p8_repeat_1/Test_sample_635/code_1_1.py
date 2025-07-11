import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and displays a left coprime factorization of the given transfer function H(s).
    """
    # Define the symbolic variable s
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # The derived left-coprime polynomial matrix D(s)
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    # Calculate the numerator polynomial matrix N(s) = D(s) * H(s)
    # sympy.simplify() is used to resolve the expressions into polynomial form
    N = sympy.simplify(D * H)

    # Display the results
    print("A left coprime factorization H(s) = D(s)^-1 * N(s) is given by:\n")
    
    print("D(s) =")
    sympy.pprint(D, use_unicode=False)
    print("\nN(s) =")
    sympy.pprint(N, use_unicode=False)
    
    # Format and print the full equation
    print("\nSo the factorization is H(s) = D(s)^-1 * N(s) where:")

    d_str = sympy.pretty(D, use_unicode=False).split('\n')
    n_str = sympy.pretty(N, use_unicode=False).split('\n')
    h_str = sympy.pretty(H, use_unicode=False).split('\n')

    # Find max widths for alignment
    max_len_d = max(len(line) for line in d_str)
    max_len_h = max(len(line) for line in h_str)

    # Print the equation with matrices aligned
    for i in range(len(d_str)):
        h_line = h_str[i] if i < len(h_str) else ''
        d_line = d_str[i] if i < len(d_str) else ''
        n_line = n_str[i] if i < len(n_str) else ''
        
        op_str = "=" if i == len(d_str) // 2 else " "
        inv_str = "^-1 *" if i == len(d_str) // 2 else "      "

        print(f"{h_line.ljust(max_len_h)} {op_str} {d_line.ljust(max_len_d)} {inv_str} {n_line}")


if __name__ == '__main__':
    calculate_left_coprime_factorization()