import math

def prime_factorization(n):
    """Computes the prime factorization of an integer n."""
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

def main():
    """
    Solves the problem by calculating the sum S and then finding the number of its divisors.
    """
    # Exponents of primes in N = 2^a1 * 29^a2 * 59^a3 * 79^a4
    a1, a2, a3, a4 = 8, 59, 79, 29

    print("Step 1: Calculate the sum components.")

    # Sum1 is for the prime 29. The term is (a2+1-b).
    # This is the sum of integers from 1 to a2+1.
    sum1 = (a2 + 1) * (a2 + 2) // 2
    print(f"The first component of the sum, Sum1, is: {sum1}")

    # Sum2 involves primes 59 and 79. We need to split based on parity of exponents.
    # Sum for even exponents of 59^c, term is (a3+1-c)
    c_even_count = a3 // 2 + 1
    c_even_sum = ((a3 + 1) + ((a3 + 1) - 2 * (c_even_count - 1))) * c_even_count // 2
    
    # Sum for odd exponents of 59^c
    c_odd_count = (a3 + 1) // 2
    c_odd_sum = (((a3 + 1) - 1) + ((a3 + 1) - (2 * (c_odd_count - 1) + 1))) * c_odd_count // 2

    # Sum for even exponents of 79^e, term is (a4+1-e)
    e_even_count = a4 // 2 + 1
    e_even_sum = ((a4 + 1) + ((a4 + 1) - 2 * (e_even_count - 1))) * e_even_count // 2

    # Sum for odd exponents of 79^e
    e_odd_count = (a4 + 1) // 2
    e_odd_sum = (((a4 + 1) - 1) + ((a4 + 1) - (2 * (e_odd_count - 1) + 1))) * e_odd_count // 2
    
    # Sum2 = (sum over even c)*(sum over even e) + (sum over odd c)*(sum over odd e)
    sum2 = c_even_sum * e_even_sum + c_odd_sum * e_odd_sum
    print(f"The second component of the sum, Sum2, is: {c_even_sum} * {e_even_sum} + {c_odd_sum} * {e_odd_sum} = {sum2}")

    # The total sum S = (a1+1) * Sum1 * Sum2
    s_val = (a1 + 1) * sum1 * sum2
    print(f"\nStep 2: Calculate the total sum S.")
    print(f"S = ({a1}+1) * {sum1} * {sum2} = {s_val}")

    print("\nStep 3: Find the prime factorization of S.")
    
    # Combine prime factorizations of each part of S
    s_factors = {}
    
    # Factors from (a1+1)
    f1 = prime_factorization(a1 + 1)
    for p, e in f1.items(): s_factors[p] = s_factors.get(p, 0) + e
    
    # Factors from Sum1
    f2 = prime_factorization(sum1)
    for p, e in f2.items(): s_factors[p] = s_factors.get(p, 0) + e

    # Factors from Sum2
    f3 = prime_factorization(sum2)
    for p, e in f3.items(): s_factors[p] = s_factors.get(p, 0) + e

    sorted_factors = sorted(s_factors.items())
    factor_str = " * ".join([f"{p}^{e}" for p, e in sorted_factors])
    print(f"The prime factorization of S is: {factor_str}")

    print("\nStep 4: Calculate the number of divisors of S.")
    num_divisors = 1
    exponents = []
    for p, e in sorted_factors:
        num_divisors *= (e + 1)
        exponents.append(e)

    calc_str = " * ".join([f"({e}+1)" for e in exponents])
    print(f"The number of divisors is the product of (exponent + 1) for each prime factor.")
    print(f"Number of divisors = {calc_str} = {num_divisors}")

if __name__ == '__main__':
    main()