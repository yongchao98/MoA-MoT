import math
from collections import Counter

def get_partitions(n):
    """Generates all integer partitions of n."""
    if n == 0:
        return [[]]
    parts = set()
    # A recursive approach to generate unique partitions
    for i in range(1, n + 1):
        for p in get_partitions(n - i):
            # Add the new part and sort to create a canonical representation for the set
            parts.add(tuple(sorted([i] + p, reverse=True)))
    # Return as a list of lists
    return [list(p) for p in parts]

def get_divisors(k):
    """Gets all divisors of an integer k."""
    divs = set()
    for i in range(1, int(math.sqrt(k)) + 1):
        if k % i == 0:
            divs.add(i)
            divs.add(k // i)
    return list(divs)

def solve():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    """
    n = 4
    n_factorial = math.factorial(n)
    
    partitions = get_partitions(n)
    # Sort partitions for a consistent calculation order, e.g., by the largest part.
    partitions.sort(key=lambda p: tuple([-x for x in p]))

    sum_terms = []

    for partition in partitions:
        # j[k] will be the number of cycles of length k in the partition.
        counts = Counter(partition)
        j = [0] * (n + 1)
        for length, num_cycles in counts.items():
            j[length] = num_cycles

        # Calculate the number of permutations with this cycle structure.
        # Formula: n! / (product_k (k^j_k * j_k!))
        perms_denominator = 1
        for k in range(1, n + 1):
            if j[k] > 0:
                perms_denominator *= (k**j[k] * math.factorial(j[k]))
        num_perms = n_factorial // perms_denominator

        # Calculate the number of functions that commute with a permutation of this structure.
        # Formula: product_k (sum_{d|k} d * j_d)^j_k
        num_commuting_funcs = 1
        for k in range(1, n + 1):
            if j[k] > 0:
                sum_term = 0
                divisors = get_divisors(k)
                for d in divisors:
                    sum_term += d * j[d]
                num_commuting_funcs *= (sum_term ** j[k])

        contribution = num_perms * num_commuting_funcs
        sum_terms.append(contribution)

    total_sum_term = sum(sum_terms)
    num_classes = total_sum_term // n_factorial

    # Print the final equation with all its components, as requested.
    final_sum_str = " + ".join(map(str, sum_terms))
    print(f"The number of classes is given by Burnside's Lemma:")
    print(f"(1 / {n_factorial}) * ({final_sum_str}) = {total_sum_term} / {n_factorial} = {num_classes}")

solve()
<<<19>>>