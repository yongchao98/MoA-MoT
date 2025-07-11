import math

def prime_factorize(n):
    """
    Returns a dictionary of prime factors of n.
    """
    factors = {}
    # Check for factor 2
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def main():
    """
    Solves the problem step-by-step.
    """
    # Exponents from N = 2^8 * 29^59 * 59^79 * 79^29
    e2 = 8
    e29 = 59
    e59 = 79
    e79 = 29

    # Part 1: Calculate the sum related to the exponent of 29
    # This sum is Sum_{a=0 to 59} (60-a)
    sum_a = sum(range(1, e29 + 2))
    
    # Part 2: Calculate sums related to exponents of 59 and 79
    # We need to sum (80-b)(30-c) over b,c where b+c is even.
    # This splits into (b even, c even) and (b odd, c odd).

    # Sums for exponent b of 59
    limit_b = e59 + 1
    sum_b_even = sum(limit_b - b for b in range(0, limit_b, 2))
    sum_b_odd = sum(limit_b - b for b in range(1, limit_b, 2))
    
    # Sums for exponent c of 79
    limit_c = e79 + 1
    sum_c_even = sum(limit_c - c for c in range(0, limit_c, 2))
    sum_c_odd = sum(limit_c - c for c in range(1, limit_c, 2))
    
    # Combine to get the sum K
    K = sum_b_even * sum_c_even + sum_b_odd * sum_c_odd

    # Part 3: Calculate the total sum S
    # S = tau(2^8) * (Sum over a) * (Sum over b,c)
    tau_2_part = e2 + 1
    S = tau_2_part * sum_a * K

    # Part 4: Find the prime factorization of S
    s_factors = prime_factorize(S)
    
    # Part 5: Calculate the number of divisors of S
    num_divisors = 1
    for p in s_factors:
        num_divisors *= (s_factors[p] + 1)
        
    # Print the results
    factor_str_parts = []
    for p in sorted(s_factors.keys()):
        factor_str_parts.append(f"{p}^{s_factors[p]}")
    print(f"The sum S has the prime factorization: {' * '.join(factor_str_parts)}")
    
    divisor_eq_parts = []
    for p in sorted(s_factors.keys()):
        divisor_eq_parts.append(f"({s_factors[p]}+1)")
    print(f"The number of divisors is {' * '.join(divisor_eq_parts)} = {num_divisors}")

if __name__ == "__main__":
    main()
