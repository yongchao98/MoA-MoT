import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator of the Jacobi symbol must be a positive odd integer.")
    
    a %= n
    if a == 0:
        return 0 if n > 1 else 1

    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                result = -result
        
        a, n = n, a
        
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
        
    if n == 1:
        return result
    else:
        return 0

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2*N1} x U(1)_{-2*N2} theory.
    'n' must be a non-zero odd integer.
    """
    if n == 0 or n % 2 == 0:
        print("zeta_n is calculated for a non-zero odd integer 'n'.")
        # For even n, the sum may be zero, making zeta_n undefined.
        # This simplified formula is valid for odd n.
        return None

    numerator = N1 * N2
    denominator = abs(n)
    
    try:
        result = jacobi_symbol(numerator, denominator)
        print(f"For n = {n}, N1 = {N1}, N2 = {N2}:")
        # Output each number in the final equation
        print(f"zeta_{n} = (({N1} * {N2}) / |{n}|) = ({numerator} / {denominator}) = {result}")
        return result
    except ValueError as e:
        print(e)
        return None

if __name__ == '__main__':
    # Example values, as if provided by the user
    N1 = 3
    N2 = 5
    n = 7 # Should be a non-zero odd integer

    # Calculate and print the result
    final_result = calculate_zeta_n(n, N1, N2)
    
    # You can also test other values
    # N1 = 4, N2 = 7, n = -3
    # calculate_zeta_n(-3, 4, 7)

<<<1>>>