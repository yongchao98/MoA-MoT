import math

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a/n).
    n must be a positive odd integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer for Jacobi symbol.")
    a %= n
    if a == 0:
        return 0
    if a == 1:
        return 1
    
    # Property: (2/n) = (-1)^((n^2-1)/8)
    if a == 2:
        if n % 8 in (1, 7):
            return 1
        else: # n % 8 in (3, 5)
            return -1

    # Property: (ab/n) = (a/n)(b/n)
    if a % 2 == 0:
        return jacobi_symbol(2, n) * jacobi_symbol(a // 2, n)

    # Law of Quadratic Reciprocity: (a/n) = (n/a) * (-1)^((a-1)(n-1)/4)
    # for odd a, n.
    if (a - 1) * (n - 1) // 4 % 2 == 0: # exponent is even
        return jacobi_symbol(n, a)
    else: # exponent is odd
        return -jacobi_symbol(n, a)

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for U(1)_{2*N1} x U(1)_{-2*N2}.
    This formula is valid for odd n.
    """
    print(f"Calculating zeta_{n} for N1={N1}, N2={N2}")
    print("-" * 30)

    if n % 2 == 0:
        print("Error: The provided formula is valid only for odd n.")
        return

    # Calculate components for the first term
    d1 = math.gcd(n, N1)
    n1 = n // d1
    m1 = N1 // d1
    
    # Calculate components for the second term
    d2 = math.gcd(n, N2)
    n2 = n // d2
    m2 = N2 // d2

    print(f"For the U(1)_{2*N1} part:")
    print(f"d1 = gcd({n}, {N1}) = {d1}")
    print(f"n1 = n/d1 = {n}/{d1} = {n1}")
    print(f"m1 = N1/d1 = {N1}/{d1} = {m1}")
    print("-" * 30)

    print(f"For the U(1)_{-2*N2} part:")
    print(f"d2 = gcd({n}, {N2}) = {d2}")
    print(f"n2 = n/d2 = {n}/{d2} = {n2}")
    print(f"m2 = N2/d2 = {N2}/{d2} = {m2}")
    print("-" * 30)

    # Calculate the phase factor
    phase_exponent = (n2 - n1) / 2
    phase_factor = 1j ** phase_exponent
    
    # Calculate the Jacobi symbols
    jacobi1 = jacobi_symbol(m1, n1)
    jacobi2 = jacobi_symbol(m2, n2)

    # Final result
    zeta_n = phase_factor * jacobi1 * jacobi2

    print("The formula is: zeta_n = i^((n2-n1)/2) * (m1/n1) * (m2/n2)")
    print(f"Plugging in the numbers:")
    print(f"zeta_{n} = i^(({n2} - {n1})/2) * ({m1}/{n1}) * ({m2}/{n2})")
    print(f"zeta_{n} = ({phase_factor}) * ({jacobi1}) * ({jacobi2})")
    print(f"zeta_{n} = {zeta_n}")

if __name__ == '__main__':
    # Example usage:
    # Please provide your own values for N1, N2, and n.
    # n must be an odd integer.
    N1_val = 3
    N2_val = 5
    n_val = 7
    calculate_zeta_n(N1_val, N2_val, n_val)
    
    print("\nAnother example:")
    N1_val = 2
    N2_val = 1
    n_val = 3
    calculate_zeta_n(N1_val, N2_val, n_val)