import math

def get_odd_part(m):
    """Returns the odd part of a positive integer."""
    if m <= 0:
        raise ValueError("Input must be a positive integer.")
    while m % 2 == 0:
        m //= 2
    return m

def jacobi_symbol(a, n):
    """
    Calculates the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator n for the Jacobi symbol must be a positive odd integer.")
    
    # Rule 2: (a/n) = (a mod n / n)
    a = a % n
    if a == 0:
        return 0
    if math.gcd(a, n) != 1:
        return 0
    
    t = 1
    while a != 0:
        # Rule 4: (2/n) part
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        # Rule 5: Law of quadratic reciprocity
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
            
        # Rule 2 again
        a = a % n
        
    if n == 1:
        return t
    else:
        # This case should not be reached if gcd(a,n)=1 initially
        return 0

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    This formula is valid for odd n, even N1, N2, and gcd(n, N1*N2)=1.
    """
    print(f"Calculating for n = {n}, N1 = {N1}, N2 = {N2}:")
    
    if n % 2 == 0:
        print("n must be an odd integer.")
        return
    if N1 % 2 != 0 or N2 % 2 != 0:
        print("This formula for zeta_n assumes N1 and N2 are even, otherwise the sum is often zero.")
        # return
        
    if math.gcd(n, N1 * N2) != 1:
        print(f"This simplified formula assumes gcd(n, N1*N2)=1. The general case is more complex.")
        # We can still compute, but mention the caveat.
    
    try:
        N1_tilde = get_odd_part(N1)
        N2_tilde = get_odd_part(N2)
        numerator = N1_tilde * N2_tilde
        
        result = jacobi_symbol(numerator, n)
        
        print(f"The odd part of N1 is N1_tilde = {N1_tilde}")
        print(f"The odd part of N2 is N2_tilde = {N2_tilde}")
        print(f"The higher central charge is given by the formula:")
        print(f"zeta_n = ( (odd_part(N1) * odd_part(N2)) / n )")
        print(f"zeta_{n} = ( ({N1_tilde} * {N2_tilde}) / {n} ) = ( {numerator} / {n} )")
        print(f"The result is: {result}")

    except ValueError as e:
        print(f"Error: {e}")


# --- User-configurable values ---
# n must be a positive odd integer
n_val = 5
# For the sum to be non-zero in the simplified case, N1 and N2 should be even.
N1_val = 12
N2_val = 22
# ------------------------------------

calculate_zeta_n(n_val, N1_val, N2_val)