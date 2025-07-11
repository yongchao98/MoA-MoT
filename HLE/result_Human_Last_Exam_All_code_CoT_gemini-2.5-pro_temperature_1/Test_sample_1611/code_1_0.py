import collections

# Grid dimensions
n = 4
m = 4

# Step 1: Find all valid column patterns for a column of height n.
# A pattern is valid if it has no adjacent 1s.
# We can check this by seeing if a bitwise AND with its shifted version is zero.
valid_patterns = []
for i in range(1 << n):
    if (i & (i >> 1)) == 0:
        valid_patterns.append(i)

# Step 2: Pre-calculate which patterns are compatible with each other.
# Two patterns p1 and p2 are compatible if placing them next to each other
# does not create adjacent 1s. This is true if their bitwise AND is 0.
compatibility_map = collections.defaultdict(list)
for p1 in valid_patterns:
    for p2 in valid_patterns:
        if (p1 & p2) == 0:
            compatibility_map[p1].append(p2)

# Step 3: Dynamic Programming to count the assignments.
# dp[p] will store the number of ways to tile an n x j grid, ending with pattern p.

# Base case (j=1): For the first column, there is 1 way for each valid pattern.
dp = {p: 1 for p in valid_patterns}

# Iterate for the remaining columns (from the 2nd to the m-th column).
for j in range(1, m):
    new_dp = collections.defaultdict(int)
    # Calculate the number of ways for each possible pattern in the current column.
    for p_new in valid_patterns:
        # This is the sum of ways for all compatible previous columns.
        count = 0
        for p_old in compatibility_map[p_new]:
            count += dp[p_old]
        new_dp[p_new] = count
    # Update dp table for the next iteration.
    dp = new_dp

# Step 4: The final answer is the sum of counts for all possible patterns in the last column.
total_assignments = sum(dp.values())

# The final result is the sum of counts for each possible configuration of the last column.
# We sort the values for a consistent output format.
final_counts = sorted(dp.values())

# Print the final breakdown of the sum as requested.
equation_parts = [str(c) for c in final_counts]
print(f"The total number of valid 0/1 assignments is the sum of counts for each possible valid last column:")
print(f"{' + '.join(equation_parts)} = {total_assignments}")

print("\nFinal Answer:")
print(f"<<<{total_assignments}>>>")