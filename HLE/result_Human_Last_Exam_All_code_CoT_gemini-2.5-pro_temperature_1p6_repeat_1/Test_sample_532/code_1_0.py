import math
from collections import defaultdict
from itertools import product

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
    factors = defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def get_partitions(k):
    """Generates all integer partitions of k."""
    if k == 0:
        yield []
        return
    for i in range(1, k + 1):
        for p in get_partitions(k - i):
            if not p or i >= p[0]:
                yield [i] + p

def get_abelian_group_structures(n):
    """
    Finds all non-isomorphic abelian groups of order n.
    This represents all finite filled nilpotent groups of order n.
    """
    if n <= 0:
        return []

    if n == 1:
        return ["{e}"]

    prime_factors = get_prime_factorization(n)
    
    # This list will hold lists of p-group structures for each prime factor
    # e.g., [[C_4, C_2xC_2], [C_3]] for n=12
    all_p_group_options = []

    sorted_primes = sorted(prime_factors.keys())

    for p in sorted_primes:
        a = prime_factors[p]
        partitions = get_partitions(a)
        
        p_group_structures = []
        for part in partitions:
            # A partition like [2, 1, 1] for p^4 means C_{p^2} x C_{p^1} x C_{p^1}
            components = [f"C_{p**val}" for val in reversed(part)]
            p_group_structures.append(" x ".join(components))
        all_p_group_options.append(p_group_structures)
        
    # Generate all combinations of p-group structures
    final_group_structures = []
    for combo in product(*all_p_group_options):
        final_group_structures.append(" x ".join(combo))
        
    return final_group_structures

if __name__ == '__main__':
    # The classification of finite filled nilpotent groups is that they are precisely the finite abelian groups.
    # This script will list all such groups for a user-specified order n.
    
    try:
        n_str = input("Enter a positive integer order n: ")
        n = int(n_str)
        if n <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print("Invalid input. Please enter a positive integer.")
    else:
        print(f"\nThe finite filled nilpotent groups are the finite abelian groups.")
        print(f"Based on this classification, here are all filled nilpotent groups of order {n}:")
        
        group_structures = get_abelian_group_structures(n)
        
        if not group_structures:
            print("No such groups found (input might be invalid or 1).")
        else:
            for group in group_structures:
                print(group)
