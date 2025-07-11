import math
from functools import reduce
from collections import defaultdict
import itertools

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = defaultdict(int)
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return dict(factors)

def get_partitions(n):
    """Generates all integer partitions of n."""
    # Using a standard recursive generator for integer partitions
    if n == 0:
        yield []
        return
    for p in get_partitions(n - 1):
        yield [1] + p
        if len(p) > 0 and p[0] > 1:
            yield [p[0] - 1] + p[1:]

def get_unique_partitions(n):
    """Returns a list of unique partitions for n."""
    return sorted([sorted(p, reverse=True) for p in get_partitions(n)])

def get_all_abelian_group_exponents(order_n):
    """
    Finds the exponents of all non-isomorphic abelian groups of a given order n.
    """
    if order_n == 1:
        return {1}

    prime_factors = get_prime_factorization(order_n)
    
    # For each prime factor p^a, find the possible exponents of the p-group component.
    # The exponent is p^max(partition) for each partition of a.
    p_group_exponent_options = []
    for p, a in prime_factors.items():
        partitions_of_a = {tuple(p) for p in get_unique_partitions(a)}
        exponents_for_p = {p ** max(part) for part in partitions_of_a}
        p_group_exponent_options.append(list(exponents_for_p))

    # The exponent of a group is the LCM of the exponents of its p-group components.
    # We find all combinations of p-group exponents.
    def lcm(a, b):
        return (a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

    group_exponents = set()
    for combo in itertools.product(*p_group_exponent_options):
        group_exponent = reduce(lcm, combo, 1)
        group_exponents.add(group_exponent)
        
    return group_exponents

def main():
    """
    Calculates the number of unique sets S(rho) cap D for representations of
    all finite Abelian groups of cardinality 18.
    """
    order = 18
    
    # As per the logic, the total number of unique sets is determined by the
    # number of divisors of the maximum possible exponent for an abelian group
    # of the given order.
    
    # 1. Find all possible exponents for abelian groups of order 18.
    group_exponents = get_all_abelian_group_exponents(order) # Should be {18, 6}
    
    # 2. Find the maximum exponent.
    max_exponent = max(group_exponents)
    
    # 3. The number of unique sets is the number of divisors of this max_exponent.
    # We calculate this using its prime factorization.
    factors = get_prime_factorization(max_exponent)
    
    # 4. Build the expression for the number of divisors and calculate it.
    # tau(n) = (a1+1)*(a2+1)*... for n = p1^a1 * p2^a2 * ...
    exponent_values = list(factors.values())
    
    # Creating the equation string to be printed
    equation_parts = [f"({exp}+1)" for exp in exponent_values]
    equation_str = " * ".join(equation_parts)
    
    # Calculating the final answer
    num_divisors = reduce(lambda x, y: x * (y + 1), exponent_values, 1)

    print(f"The number of unique sets is the number of divisors of {max_exponent}.")
    print("This is calculated from its prime factorization:")
    prime_factor_str = " * ".join([f"{p}^{a}" for p, a in factors.items()])
    print(f"{max_exponent} = {prime_factor_str}")
    print("The number of divisors is given by the equation:")
    print(f"{equation_str} = {num_divisors}")


if __name__ == "__main__":
    main()
>>>6