import math
from collections import Counter

def combinations_with_replacement(n, k):
    """
    Calculates the number of ways to choose k items from n types with replacement.
    This is equivalent to C(n + k - 1, k).
    """
    if k < 0 or n < 0:
        return 0
    if k == 0 or n == 0:
        return 1 if n > 0 or k == 0 else 0
    return math.comb(n + k - 1, k)

def generate_partitions(n):
    """Generates all unique integer partitions of n."""
    partitions_set = set()
    def find_partitions(target, current_partition):
        if target == 0:
            partitions_set.add(tuple(sorted(current_partition, reverse=True)))
            return
        
        start = current_partition[-1] if current_partition else 1
        
        for i in range(start, target + 1):
            find_partitions(target - i, current_partition + [i])

    find_partitions(n, [])
    return sorted(list(partitions_set), key=lambda p: (len(p), p), reverse=True)

def solve():
    """
    Calculates the number of non-isomorphic endomorphisms on a set of size 4.
    """
    n = 4
    
    # C(k): Number of connected functional graphs on k vertices.
    C = {
        1: 1,  # 1-cycle
        2: 2,  # 2-cycle OR rooted tree on 2
        3: 4,  # 3-cycle OR 2-cycle w/ tail OR rooted tree on 3 (2 types)
        4: 9   # 4-cycle OR 3-cycle w/ tail OR 2-cycle w/ 2 tails (3 types) OR rooted tree on 4 (4 types)
    }
    
    print("The problem is equivalent to counting the number of non-isomorphic functional graphs on 4 vertices.")
    print("This is found by summing the possibilities over all integer partitions of 4.")
    print("\nFirst, the number of connected functional graphs C(k) for k=1, 2, 3, 4 are:")
    print("C(1) = 1")
    print("C(2) = 2")
    print("C(3) = 4")
    print("C(4) = 9")
    
    partitions_of_4 = generate_partitions(n)
    total_count = 0
    
    print("\nCalculating the number of structures for each partition of 4:")
    
    results_per_partition = []
    for p in partitions_of_4:
        counts = Counter(p)
        term_result = 1
        
        for size, num_occurrences in counts.items():
            num_choices = combinations_with_replacement(C[size], num_occurrences)
            term_result *= num_choices
        
        results_per_partition.append(term_result)
        total_count += term_result
        
        partition_str = "+".join(map(str, p))
        print(f"Partition {partition_str}: {term_result} structures")

    print("\nThe final calculation is the sum of the results for each partition:")
    sum_str = " + ".join(map(str, results_per_partition))
    print(f"{sum_str} = {total_count}")

solve()
<<<19>>>