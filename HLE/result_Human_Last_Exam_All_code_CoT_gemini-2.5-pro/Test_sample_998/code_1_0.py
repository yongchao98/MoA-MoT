import math
from collections import Counter

def prime_factorize(n):
    """
    Returns a Counter object with the prime factorization of n.
    """
    factors = Counter()
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def main():
    """
    Solves the number theory problem step-by-step.
    """
    # Step 1: Define the exponents from N = 2^8 * 29^59 * 59^79 * 79^29
    A = 8
    a1 = 59
    c1 = 79
    c2 = 29

    # Step 2: Calculate the three components of the sum S(N)
    # Component from the power of 2
    s_2 = A + 1
    
    # Component from primes = 1 (mod 4), i.e., 29^59
    s_p = (a1 + 1) * (a1 + 2) // 2
    
    # Component from primes = 3 (mod 4), i.e., 59^79 * 79^29
    # T1 and T2 are sums over all exponents
    T1 = (c1 + 1) * (c1 + 2) // 2
    T2 = (c2 + 1) * (c2 + 2) // 2
    
    # D1 and D2 are alternating sums over exponents
    # As c1 and c2 are odd, the formula is (c+1)/2
    D1 = (c1 + 1) // 2
    D2 = (c2 + 1) // 2
    
    s_q = (T1 * T2 + D1 * D2) // 2

    # Step 3: Calculate the total sum S(N)
    S_N_val = s_2 * s_p * s_q
    
    # Step 4: Find the prime factorization of S(N) by factorizing its components
    factors_s2 = prime_factorize(s_2)
    factors_sp = prime_factorize(s_p)
    factors_sq = prime_factorize(s_q)
    
    total_factors = factors_s2 + factors_sp + factors_sq

    # Step 5: Calculate the number of divisors of S(N)
    num_divisors = 1
    equation_terms = []
    
    # Sort primes for consistent output order
    sorted_primes = sorted(total_factors.keys())

    for prime in sorted_primes:
        exponent = total_factors[prime]
        num_divisors *= (exponent + 1)
        equation_terms.append(f"({exponent}+1)")

    # Print the final calculation as requested
    print("The number of divisors is the product of (exponent + 1) for each prime factor.")
    print(f"The prime factorization of the sum has exponents: {list(total_factors[p] for p in sorted_primes)}")
    print("Final calculation:")
    print(f"{' * '.join(equation_terms)} = {num_divisors}")


if __name__ == "__main__":
    main()
