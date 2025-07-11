import math

def solve_weyl_group_count():
    """
    Calculates the number of non-isomorphic finite Weyl groups of rank 4.

    This is done by considering all possible partitions of the rank 4 into
    ranks of irreducible components.
    """
    print("This script calculates the number of non-isomorphic finite Weyl groups of rank 4.")
    print("The calculation is based on the classification of irreducible finite Weyl groups and partitions of the rank.\n")

    # The classification of irreducible finite Weyl groups by rank.
    # Note: Weyl groups for B_n and C_n are isomorphic, so we only list B_n.
    irreducible_by_rank = {
        1: ["A_1"],
        2: ["A_2", "B_2", "G_2"],
        3: ["A_3", "B_3"],
        4: ["A_4", "B_4", "D_4", "F_4"]
    }

    # We will sum the counts for each partition of 4.
    counts = []

    # Case 1: Irreducible groups (Partition: 4)
    # These are the connected Dynkin diagrams with 4 nodes.
    count_4 = len(irreducible_by_rank[4])
    counts.append(count_4)
    print(f"Number of irreducible groups of rank 4: {count_4}")

    # Case 2: Reducible groups (Partition: 3 + 1)
    # A product of one rank-3 group and one rank-1 group.
    count_3_1 = len(irreducible_by_rank[3]) * len(irreducible_by_rank[1])
    counts.append(count_3_1)
    print(f"Number of reducible groups for partition 3+1: {count_3_1}")

    # Case 3: Reducible groups (Partition: 2 + 2)
    # A product of two rank-2 groups. This is a combination with repetition.
    # We choose 2 groups from the set of rank-2 groups with size n=3.
    # The formula is C(n+k-1, k) where n=3, k=2. C(3+2-1, 2) = C(4, 2) = 6.
    n = len(irreducible_by_rank[2])
    k = 2
    count_2_2 = math.comb(n + k - 1, k)
    counts.append(count_2_2)
    print(f"Number of reducible groups for partition 2+2: {count_2_2}")

    # Case 4: Reducible groups (Partition: 2 + 1 + 1)
    # A product of one rank-2 group and two rank-1 groups.
    # Since there is only one type of rank-1 group (A_1), the number of
    # combinations is just the number of choices for the rank-2 group.
    count_2_1_1 = len(irreducible_by_rank[2])
    counts.append(count_2_1_1)
    print(f"Number of reducible groups for partition 2+1+1: {count_2_1_1}")

    # Case 5: Reducible groups (Partition: 1 + 1 + 1 + 1)
    # A product of four rank-1 groups. Since there's only one type of rank-1
    # group, there is only one such combination: A_1 x A_1 x A_1 x A_1.
    count_1_1_1_1 = 1
    counts.append(count_1_1_1_1)
    print(f"Number of reducible groups for partition 1+1+1+1: {count_1_1_1_1}")

    # Final summation
    total_count = sum(counts)
    equation = " + ".join(map(str, counts))
    print("\nTotal number is the sum of these cases:")
    print(f"Final Equation: {equation} = {total_count}")

solve_weyl_group_count()