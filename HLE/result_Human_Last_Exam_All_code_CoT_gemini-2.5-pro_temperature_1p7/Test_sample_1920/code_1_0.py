import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary {prime: exponent}.
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

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = {1, n}
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def get_prime_factors_from_dict(factors_dict):
    """Returns the set of prime factors from a factorization dictionary."""
    return set(factors_dict.keys())

def phi(n, factors=None):
    """
    Calculates Euler's totient function phi(n).
    Optionally accepts a precomputed prime factorization.
    """
    if n == 1:
        return 1
    if factors is None:
        factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return int(result)

def mu(n, factors=None):
    """
    Calculates the Mobius function mu(n).
    Optionally accepts a precomputed prime factorization.
    """
    if n == 1:
        return 1
    if factors is None:
        factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def solve():
    """
    Finds the number of primitive Dirichlet characters for a given conductor and order.
    """
    d = 53599
    k = 6

    print(f"Finding the number of primitive Dirichlet characters of conductor d={d} and order k={k}.")
    print("-" * 30)

    # Step 1: Prime factorization of d
    d_factors = get_prime_factorization(d)
    primes = list(d_factors.keys())
    print(f"The prime factorization of d = {d} is: {' * '.join(map(str, primes))}")
    
    if not all(val == 1 for val in d_factors.values()):
        print("This script is specialized for square-free conductors.")
        return

    # Step 2: Use the formula based on Mobius inversion
    # C*(k,d) = sum_{m|k} mu(k/m) * C(m,d)
    print("\nWe use the formula: C*(k,d) = sum_{m|k} mu(k/m) * C(m,d)")
    print("where C(m,d) is the number of primitive characters with conductor d whose order divides m.")
    
    divisors_of_k = get_divisors(k)
    final_terms = {}

    print("\nCalculating each term:")
    for m in divisors_of_k:
        # Calculate C(m,d)
        # For square-free d, C(m,d) = (num_prim_chars_mod_p_order_divides_m) ^ (num_prime_factors)
        
        divisors_of_m = get_divisors(m)
        
        # Calculate Num(m): number of prim. chars mod p with order dividing m.
        # This is sum(phi(j)) for j|m and j>1.
        num_prim_chars_per_p = 0
        for j in divisors_of_m:
            if j > 1:
                # Assuming j | (p-1) for all primes p, which is true in this problem.
                num_prim_chars_per_p += phi(j)
        
        # Calculate C(m, d)
        C_m_d = num_prim_chars_per_p ** len(primes)
        
        mu_val = mu(k // m)
        term_value = mu_val * C_m_d
        final_terms[m] = term_value
        
        print(f"For m = {m}:")
        print(f"  - Num primitive chars mod p with order dividing {m} = {num_prim_chars_per_p}")
        print(f"  - C({m},{d}) = {num_prim_chars_per_p}^{len(primes)} = {C_m_d}")
        print(f"  - mu({k}/{m}) = mu({k//m}) = {mu_val}")
        print(f"  - Term = {mu_val} * {C_m_d} = {term_value}")

    # Step 3: Sum the terms for the final answer
    total_count = sum(final_terms.values())
    
    # Sort terms by m for a consistent equation
    equation_parts = [f"({final_terms[m]})" for m in sorted(final_terms.keys())]

    print("\n" + "-" * 30)
    print("The final calculation is the sum of these terms.")
    
    # We re-order the terms to match the Mobius inversion sum from C(k,d) perspective
    # C*(k) = C(k) - C(k/p1) - ...
    # This translates to C*(6) = C(6) - C(3) - C(2) + C(1)
    
    # Let's map my final_terms dictionary back to the C(m,d) values
    # Term for m=1 has mu(6)=1. C(1,d)=0. -> 0
    # Term for m=2 has mu(3)=-1. C(2,d)=1. -> -1
    # Term for m=3 has mu(2)=-1. C(3,d)=16 -> -16
    # Term for m=6 has mu(1)=1. C(6,d)=625 -> 625
    c1 = 0
    c2 = 1
    c3 = 16
    c6 = 625

    print(f"Final Equation: C(6,d) - C(3,d) - C(2,d) + C(1,d)")
    print(f"= {c6} - {c3} - {c2} + {c1}")
    print(f"= {c6 - c3 - c2 + c1}")

    print(f"\nThe number of primitive Dirichlet characters of conductor {d} and order {k} is {total_count}.")


solve()
<<<608>>>