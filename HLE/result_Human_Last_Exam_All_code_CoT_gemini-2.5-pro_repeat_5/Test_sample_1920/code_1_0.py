import math
from itertools import combinations

def prime_factorize(n):
    """Returns a dictionary of prime factors of n."""
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
    """Returns a set of divisors of n."""
    divs = {1}
    i = 2
    while i * i <= n:
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
        i += 1
    if n > 1:
        divs.add(n)
    return divs

_phi_cache = {}
def euler_phi(n):
    """Computes Euler's totient function."""
    if n in _phi_cache:
        return _phi_cache[n]
    if n == 1:
        _phi_cache[1] = 1
        return 1
    factors = prime_factorize(n)
    result = n
    for p in factors:
        result -= result // p
    _phi_cache[n] = result
    return result

def solve():
    """
    Calculates the number of primitive Dirichlet characters of a given conductor and order.
    """
    d = 53599
    g = 6

    d_factors = prime_factorize(d)
    g_factors = prime_factorize(g)

    # Check for square-free conductor
    if any(a > 1 for a in d_factors.values()):
        print("The conductor d is not square-free. This script handles square-free conductors.")
        return

    d_primes = sorted(list(d_factors.keys()))
    g_primes = sorted(list(g_factors.keys()))

    print(f"Finding the number of primitive Dirichlet characters of conductor d = {d} and order g = {g}.")
    print(f"The conductor d is square-free with prime factors: {d_primes}")
    print(f"The prime factors of the order g = {g} are: {g_primes}")

    g_divisors = sorted([div for div in get_divisors(g) if div > 1])

    # Using the Principle of Inclusion-Exclusion
    terms = {}
    
    # Iterate through the power set of g_primes
    num_g_primes = len(g_primes)
    for i in range(1 << num_g_primes):
        subset_g_primes = []
        for j in range(num_g_primes):
            if (i >> j) & 1:
                subset_g_primes.append(g_primes[j])
        
        # Calculate the number of character tuples where orders are NOT divisible by any prime in subset_g_primes.
        term_count = 1
        for p in d_primes:
            count_per_prime = 0
            for k in g_divisors:
                # Check if k's prime factors are disjoint from subset_g_primes
                is_disjoint = True
                for q in subset_g_primes:
                    if k % q == 0:
                        is_disjoint = False
                        break
                
                if is_disjoint and (p - 1) % k == 0:
                    count_per_prime += euler_phi(k)
            term_count *= count_per_prime
        
        terms[tuple(sorted(subset_g_primes))] = term_count
    
    # Calculate the final result using I-E
    total = terms[()]
    subtractions = [terms[(q,)] for q in g_primes]
    additions_len_2 = [terms[tuple(c)] for c in combinations(g_primes, 2)]
    
    final_result = total - sum(subtractions) + sum(additions_len_2)

    print("\nUsing the Principle of Inclusion-Exclusion:")
    print(f"1. Total tuples of non-trivial characters with order dividing 6: {total}")
    print(f"2. Subtract tuples where order is not divisible by 2 (i.e., order is 3): {subtractions[0]}")
    print(f"3. Subtract tuples where order is not divisible by 3 (i.e., order is 2): {subtractions[1]}")
    if additions_len_2 and additions_len_2[0] > 0:
        print(f"4. Add back tuples where order is not divisible by 2 or 3 (i.e., order is 1): {additions_len_2[0]}")
        equation_str = f"{total} - {subtractions[0]} - {subtractions[1]} + {additions_len_2[0]} = {final_result}"
    else:
        # Case where this intersection is empty, which it is here since we require order > 1.
        equation_str = f"{total} - {subtractions[0]} - {subtractions[1]} + 0 = {final_result}"

    print(f"\nThe number of primitive Dirichlet characters is given by the equation:")
    print(equation_str)

solve()