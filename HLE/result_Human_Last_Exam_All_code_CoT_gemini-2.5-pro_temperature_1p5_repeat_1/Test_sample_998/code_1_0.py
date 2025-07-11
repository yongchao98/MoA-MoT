import math
from collections import defaultdict

def calculate_sum(limit, val_offset):
    """
    Calculates three sums for an exponent with range [0, limit]:
    1. The total sum of (val_offset - exponent).
    2. The sum for even-indexed exponents.
    3. The sum for odd-indexed exponents.
    """
    total_sum = (limit + 1) * (val_offset + val_offset - limit) // 2
    
    # Calculate sums for even and odd indices
    # This is equivalent to sum_{k=0..limit} (val_offset-k)x^k for x=1 and x=-1
    if (limit + 1) % 2 == 0: # even number of terms
        alternating_sum = (limit + 1) // 2
    else: # odd number of terms
        alternating_sum = val_offset - limit

    sum_even = (total_sum + alternating_sum) // 2
    sum_odd = (total_sum - alternating_sum) // 2
    return total_sum, sum_even, sum_odd

def get_prime_factorization(n, factors):
    """
    Finds the prime factorization of n and adds it to the factors dictionary.
    """
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1

def main():
    """
    Solves the problem by calculating S and then its number of divisors.
    """
    # Exponents of N = 2^8 * 29^59 * 59^79 * 79^29
    a0 = 8
    a1_lim = 59
    a2_lim = 79
    a3_lim = 29

    # Calculate sum over b (exponent of 29)
    sum_b, _, _ = calculate_sum(a1_lim, a1_lim + 1)
    
    # Calculate sums over c (exponent of 59)
    _, sum_c_even, sum_c_odd = calculate_sum(a2_lim, a2_lim + 1)

    # Calculate sums over e (exponent of 79)
    _, sum_e_even, sum_e_odd = calculate_sum(a3_lim, a3_lim + 1)

    # Sum over c and e, where c+e is even
    sum_ce = (sum_c_even * sum_e_even) + (sum_c_odd * sum_e_odd)

    # Total sum S
    term1 = a0 + 1
    term2 = sum_b
    term3 = sum_ce
    
    # We find the prime factorization of S by factoring its three integer components
    s_factors = defaultdict(int)
    get_prime_factorization(term1, s_factors)
    get_prime_factorization(term2, s_factors)
    get_prime_factorization(term3, s_factors)
    
    # Calculate the number of divisors from the prime factorization
    num_divisors = 1
    
    # Constructing the output strings
    factor_strs = []
    final_eq_nums = []
    
    sorted_primes = sorted(s_factors.keys())
    
    for p in sorted_primes:
        exponent = s_factors[p]
        factor_strs.append(f"{p}^{exponent}")
        final_eq_nums.append(f"({exponent}+1)")
        num_divisors *= (exponent + 1)

    print(f"The sum S is the product of three terms: {term1}, {term2}, and {term3}.")
    print(f"The prime factorization of S is {' * '.join(factor_strs)}.")
    print(f"The number of divisors of S is {'*'.join(final_eq_nums)} = {num_divisors}")

if __name__ == "__main__":
    main()