import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator n for the Jacobi symbol must be a positive odd integer.")
    
    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        a, n = n, a
        
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
        
    if n == 1:
        return t
    else:
        return 0

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2N1} x U(1)_{-2N2} theory.
    This formula assumes N1 and N2 are odd integers.
    """
    try:
        # The denominator for the Jacobi symbol must be odd.
        denominator = N1 * N2
        if denominator % 2 == 0:
             print(f"Error: The formula requires N1 and N2 to be odd, but N1*N2 = {denominator} is even.")
             return None

        # The derived formula for zeta_n is the Jacobi symbol (n / (N1*N2)).
        result = jacobi_symbol(n, denominator)
        return result

    except ValueError as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    # Example values, you can change them.
    # We assume n is coprime to N1*N2 as is common in this context,
    # though the Jacobi symbol function can handle the case where they share factors.
    n_val = 3
    N1_val = 5
    N2_val = 7

    print(f"Calculating the higher central charge zeta_n for the theory U(1)_(2*N1) x U(1)_(-2*N2)")
    print(f"Using the parameters:")
    print(f"n = {n_val}")
    print(f"N1 = {N1_val}")
    print(f"N2 = {N2_val}")
    
    zeta_n_value = calculate_zeta_n(n_val, N1_val, N2_val)
    
    if zeta_n_value is not None:
        # The final equation is zeta_n = (n / (N1*N2))
        print("\nThe final equation is zeta_n = (n / (N1 * N2))")
        print(f"zeta_{n_val} = ({n_val} / ({N1_val} * {N2_val})) = ({n_val} / {N1_val * N2_val})")
        print(f"The result is: {zeta_n_value}")
