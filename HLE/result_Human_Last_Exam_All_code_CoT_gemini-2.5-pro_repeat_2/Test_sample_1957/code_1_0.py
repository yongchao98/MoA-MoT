# The degree of the polynomial defining the variety S is d=6.
d = 6

# We are looking for positive integers m and n such that m*n = d.
# We want to minimize the expression 1000*m + n.

# We list all possible pairs of (m, n)
possible_pairs = []
for m in range(1, d + 1):
    if d % m == 0:
        n = d // m
        possible_pairs.append((m, n))

# Calculate the value for each pair and find the minimum
min_val = float('inf')
best_pair = (0, 0)

print("Possible (m, n) pairs and the corresponding values of 1000m + n:")
for m, n in possible_pairs:
    val = 1000 * m + n
    print(f"For m = {m}, n = {n}: 1000 * {m} + {n} = {val}")
    if val < min_val:
        min_val = val
        best_pair = (m, n)

print(f"\nThe minimum value is {min_val}, which occurs for m={best_pair[0]} and n={best_pair[1]}.")