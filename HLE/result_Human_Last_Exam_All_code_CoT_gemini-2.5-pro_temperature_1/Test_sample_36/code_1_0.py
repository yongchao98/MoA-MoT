import math
from functools import reduce
from itertools import product

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
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

def get_partitions(n):
    """Generates all integer partitions of n."""
    # Using a standard recursive generator for partitions
    if n == 0:
        yield []
        return
    for p in get_partitions(n - 1):
        yield [1] + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield [p[0] + 1] + p[1:]

def get_exponents_of_abelian_groups(n):
    """
    Calculates the exponents of all non-isomorphic abelian groups of order n.
    """
    if n == 1:
        return [1]
        
    prime_factors = get_prime_factorization(n)
    
    # For each prime factor p^a, find partitions of the exponent 'a'.
    # Each partition corresponds to a way to form p-groups.
    # e.g., for 3^2, partitions of 2 are [2] -> C_9 and [1,1] -> C_3 x C_3
    per_prime_partitions = []
    for p, a in prime_factors.items():
        p_groups = []
        # Use a set to avoid duplicates from the generator
        partitions_set = {tuple(sorted(p)) for p in get_partitions(a)}
        for part in partitions_set:
            p_groups.append([p**i for i in part])
        per_prime_partitions.append(p_groups)

    # Combine the p-group structures for different primes
    # using the Chinese Remainder Theorem logic.
    group_factor_combos = product(*per_prime_partitions)
    
    exponents = []
    
    def lcm(a, b):
        return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

    for combo in group_factor_combos:
        all_factors = []
        for factor_list in combo:
            all_factors.extend(factor_list)
        
        # The exponent of the group is the LCM of the orders of its cyclic factors.
        group_exponent = reduce(lcm, all_factors, 1)
        exponents.append(group_exponent)
        
    return exponents

def count_divisors(n):
    """Counts the number of divisors of n (tau function)."""
    if n == 0:
        return 0
    if n == 1:
        return 1, {1:1}, ["(0 + 1)"], ["1"]
        
    prime_factors = get_prime_factorization(n)
    count = 1
    
    factor_strings = []
    term_strings = []
    val_strings = []

    for p, a in sorted(prime_factors.items()):
        count *= (a + 1)
        factor_strings.append(f"{p}^{a}")
        term_strings.append(f"({a} + 1)")
        val_strings.append(str(a+1))
        
    return count, prime_factors, term_strings, val_strings

def solve():
    """
    Main function to solve the problem.
    """
    order = 18
    
    # 1. Get exponents of all non-isomorphic abelian groups of order 18
    exponents = get_exponents_of_abelian_groups(order)
    print(f"The exponents of the non-isomorphic Abelian groups of order {order} are: {sorted(exponents)}")
    
    # 2. Find the LCM of these exponents
    def lcm(a, b):
        return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

    master_exponent = reduce(lcm, exponents, 1)
    print(f"The least common multiple (LCM) of these exponents is: {master_exponent}")
    
    # 3. Count the divisors of the master exponent
    num_divisors, factors, term_strs, val_strs = count_divisors(master_exponent)
    
    factorization_str = " * ".join([f"{p}^{a}" for p, a in sorted(factors.items())])
    print(f"The prime factorization of {master_exponent} is: {factorization_str}")

    # 4. Print the final calculation and answer
    calc_str = " * ".join(term_strs)
    val_str = " * ".join(val_strs)
    
    print("The number of unique sets is the number of divisors of the LCM.")
    print(f"This is calculated as: {calc_str} = {val_str} = {num_divisors}")
    
    print(f"\n<<< {num_divisors} >>>")

solve()