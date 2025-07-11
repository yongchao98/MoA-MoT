import math

def solve_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    
    The calculation is based on the classification of finite Weyl groups.
    A finite Weyl group is a direct product of irreducible Weyl groups, and
    its rank is the sum of the ranks of its factors. We need to find all
    combinations of irreducible groups whose ranks sum to 4.

    The irreducible Weyl groups are classified by rank:
    - Rank 1: A1
    - Rank 2: A2, B2, G2
    - Rank 3: A3, B3
    - Rank 4: A4, B4, D4, F4
    """

    # Lists of irreducible Weyl groups by rank
    rank1_types = ["A1"]
    rank2_types = ["A2", "B2", "G2"]
    rank3_types = ["A3", "B3"]
    rank4_types = ["A4", "B4", "D4", "F4"]

    # We calculate the number of combinations for each partition of 4.

    # Partition 1: A single group of rank 4
    # The groups are A4, B4, D4, F4.
    count_4 = len(rank4_types)

    # Partition 2: A group of rank 3 and a group of rank 1 (3+1)
    # Combinations are (A3 x A1), (B3 x A1).
    count_3_1 = len(rank3_types) * len(rank1_types)

    # Partition 3: Two groups of rank 2 (2+2)
    # This is a combination with replacement problem. We choose 2 groups
    # from the set of rank 2 types {A2, B2, G2}.
    # The formula for combinations with replacement is C(n+k-1, k).
    n = len(rank2_types)
    k = 2
    count_2_2 = math.comb(n + k - 1, k)

    # Partition 4: One group of rank 2 and two groups of rank 1 (2+1+1)
    # The two rank 1 groups must be A1 x A1. We just need to choose one rank 2 group.
    count_2_1_1 = len(rank2_types)

    # Partition 5: Four groups of rank 1 (1+1+1+1)
    # There is only one possibility: A1 x A1 x A1 x A1.
    count_1_1_1_1 = 1
    
    # Calculate the total number of non-isomorphic groups
    total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

    # Print the explanation and the final equation
    print("The number of non-isomorphic finite Weyl groups of rank 4 is the sum of the counts for each partition of 4:")
    print(f"Partition (4): {count_4} groups")
    print(f"Partition (3+1): {count_3_1} groups")
    print(f"Partition (2+2): {count_2_2} groups")
    print(f"Partition (2+1+1): {count_2_1_1} groups")
    print(f"Partition (1+1+1+1): {count_1_1_1_1} group")
    print("\nTotal count as a sum:")
    print(f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")

solve_weyl_groups_rank_4()