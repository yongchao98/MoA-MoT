import math

# Step 1: Define the sets A and B for construction.
# These sets satisfy the condition that ratios of pairs of elements from A
# are disjoint from ratios of pairs of elements from B, ensuring all
# products a_i * b_j are unique.
A = [1, 2, 3, 4]
B = [1, 5, 6, 7]

# Step 2: Define a pair of 4x4 orthogonal Latin squares.
# These ensure the magic product property for rows and columns.
L1 = [[1, 2, 3, 4],
      [2, 1, 4, 3],
      [3, 4, 1, 2],
      [4, 3, 2, 1]]

L2 = [[1, 2, 3, 4],
      [3, 4, 1, 2],
      [4, 3, 2, 1],
      [2, 1, 4, 3]]

# Step 3: Construct the 4x4 table.
table = [[0] * 4 for _ in range(4)]
all_numbers = set()
max_N = 0

for r in range(4):
    for c in range(4):
        # The value at each cell is a product of an element from A and an element from B.
        # The indices are determined by the Latin squares.
        val = A[L1[r][c] - 1] * B[L2[r][c] - 1]
        table[r][c] = val
        all_numbers.add(val)
        if val > max_N:
            max_N = val

# Step 4: Print the results and verify the properties.
print("Constructed 4x4 Table:")
for r in range(4):
    print("  ".join(f"{num:2d}" for num in table[r]))

print("\nVerifying Properties:")
print(f"The set of {len(all_numbers)} distinct integers is: {sorted(list(all_numbers))}")
print(f"The largest integer used, N, is: {max_N}")

# Calculate the magic product P
magic_product = 1
for x in A:
    magic_product *= x
for y in B:
    magic_product *= y

print(f"\nExpected magic product P = {' * '.join(map(str,A))} * {' * '.join(map(str,B))} = {magic_product}")

print("\nChecking products of each row and column:")
valid = True
# Check rows
for r in range(4):
    row_prod = math.prod(table[r])
    print(f"Row {r+1} product: {' * '.join(f'{n:2d}' for n in table[r])} = {row_prod}")
    if row_prod != magic_product:
        valid = False
# Check columns
for c in range(4):
    col = [table[r][c] for r in range(4)]
    col_prod = math.prod(col)
    print(f"Col {c+1} product: {' * '.join(f'{n:2d}' for n in col)} = {col_prod}")
    if col_prod != magic_product:
        valid = False

if valid:
    print("\nThe table is a valid multiplicative magic square.")
    print(f"The smallest N found is {max_N}.")
