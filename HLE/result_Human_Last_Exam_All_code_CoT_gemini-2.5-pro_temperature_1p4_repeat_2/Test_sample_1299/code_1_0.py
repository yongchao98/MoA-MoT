import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors and their exponents.
    e.g., get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def sigma_0(n):
    """
    Calculates the number of divisors of n, also known as tau(n).
    """
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    num_divisors = 1
    for p in factors:
        num_divisors *= (factors[p] + 1)
    return num_divisors

def calculate_cardinalities(ell):
    """
    Calculates |U_l| and |T_l| for a given positive integer l.
    """
    if not isinstance(ell, int) or ell <= 0:
        print("Please provide a positive integer for l.")
        return

    print(f"Calculating for l = {ell}")
    
    # Get prime factorization for l
    factors_l = get_prime_factorization(ell)
    
    # --- Calculate |U_l| ---
    # Formula: |U_l| = (v2(l) + 1) * sigma_0((l / 2**v2(l))**2)
    v2_l = factors_l.get(2, 0)
    m = ell // (2**v2_l)
    m_squared = m**2
    sigma_0_m_squared = sigma_0(m_squared)
    
    U_l = (v2_l + 1) * sigma_0_m_squared
    
    print("\n--- Part A: Calculating |U_l| ---")
    print(f"Prime factorization of {ell}: {factors_l}")
    print(f"Exponent of 2 in l (v2(l)): {v2_l}")
    print(f"Odd part of l (m): {m}")
    print(f"Number of divisors of m^2, sigma_0({m_squared}): {sigma_0_m_squared}")
    print(f"|U_{ell}| = (v2(l) + 1) * sigma_0(m^2) = ({v2_l} + 1) * {sigma_0_m_squared} = {U_l}")

    # --- Calculate |T_l| ---
    # Formula: |T_l| = sigma_0(l^2) - 1 for l > 1, and 1 for l=1
    if ell == 1:
        T_l = 1
        print("\n--- Part B (Original): Calculating |T_l| ---")
        print(f"|T_1| = 1")

    else:
        l_squared = ell**2
        sigma_0_l_squared = sigma_0(l_squared)
        T_l = sigma_0_l_squared - 1
        
        print("\n--- Part B (Original): Calculating |T_l| ---")
        print(f"l^2 = {l_squared}")
        print(f"Number of divisors of l^2, sigma_0({l_squared}): {sigma_0_l_squared}")
        print(f"|T_{ell}| = sigma_0(l^2) - 1 = {sigma_0_l_squared} - 1 = {T_l}")

# Example calculation for l = 12
calculate_cardinalities(12)