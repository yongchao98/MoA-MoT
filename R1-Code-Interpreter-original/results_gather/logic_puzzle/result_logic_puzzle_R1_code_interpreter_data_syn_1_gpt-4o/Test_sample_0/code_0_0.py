# Define the grid
grid = [
    [1, 2],
    [3, 1]
]

# Selected positions
selected_positions = [(0, 1), (0, 0), (1, 1)]

# Calculate row and column sums
row_sums = [0, 0]
col_sums = [0, 0]

for r, c in selected_positions:
    row_sums[r] += grid[r][c]
    col_sums[c] += grid[r][c]

# Check if all sums are <= 4
valid = all(sum_val <= 4 for sum_val in row_sums + col_sums)

# Output the result
print(selected_positions if valid else "No valid selection")