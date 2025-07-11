import math
from collections import Counter

def get_prime_factorization(n):
    """
    Returns a Counter object representing the prime factorization of n.
    """
    factors = Counter()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return factors

def main():
    """
    This function solves the problem step-by-step as outlined above.
    """
    # Let N = 2^8 * 29^59 * 59^79 * 79^29.
    # We want to find the number of divisors of S = sum_{d|N} f(d).
    # S = sum_{c|N, c=1(mod 4)} tau(N/c)
    # S = 9 * (sum for e1) * (sum for e2, e3)

    print("Step 1: Calculate the components of the sum S.")
    
    # This is the (8+1) part from tau(N/c) as e0=0.
    c0 = 9
    
    # Calculate Sum 1: sum_{e1=0 to 59} (60-e1) for prime 29
    # The exponent range for 29 is 0 to 59.
    s1 = sum(60 - e1 for e1 in range(60))

    # Calculate Sum 2: This involves primes 59 and 79.
    # Let A be the sum part for prime 59 (exponent e2), B for prime 79 (exponent e3)
    # Exponent range for 59 is 0 to 79.
    # Exponent range for 79 is 0 to 29.
    
    # Sum over even exponents for 59
    sum_A_even = sum(80 - e2 for e2 in range(0, 80, 2))
    # Sum over odd exponents for 59
    sum_A_odd = sum(80 - e2 for e2 in range(1, 80, 2))
    
    # Sum over even exponents for 79
    sum_B_even = sum(30 - e3 for e3 in range(0, 30, 2))
    # Sum over odd exponents for 79
    sum_B_odd = sum(30 - e3 for e3 in range(1, 30, 2))
    
    # The sum T is for pairs (e2, e3) with same parity
    # T = (e2 even, e3 even) + (e2 odd, e3 odd)
    s2 = sum_A_even * sum_B_even + sum_A_odd * sum_B_odd
    
    # The total sum S
    S = c0 * s1 * s2
    
    print(f"The total sum S is an equation of 3 parts: C0 * S1 * S2")
    print(f"Part C0 (from prime 2) = {c0}")
    print(f"Part S1 (from prime 29) = {s1}")
    print(f"Part S2 (from primes 59, 79) = {s2}")
    print(f"S = {c0} * {s1} * {s2} = {S}\n")
    
    print("Step 2: Find the prime factorization of S.")
    total_factors = get_prime_factorization(S)
    sorted_factors = sorted(total_factors.items())
    
    factor_str = ' * '.join([f'{p}^{e}' for p, e in sorted_factors])
    print(f"The prime factorization of S is: {factor_str}\n")
    
    print("Step 3: Calculate the number of divisors from the factorization.")
    num_divisors = 1
    for p, e in sorted_factors:
        num_divisors *= (e + 1)
        
    calc_str = ' * '.join([f'({e}+1)' for p, e in sorted_factors])
    vals_str = ' * '.join([str(e + 1) for p, e in sorted_factors])
    
    print(f"The number of divisors is the product of (exponent + 1) for each prime factor.")
    print(f"Number of divisors = {calc_str} = {vals_str} = {num_divisors}")
    
    print(f"\nThe final answer is {num_divisors}.")
    print(f"<<<{num_divisors}>>>")

if __name__ == "__main__":
    main()