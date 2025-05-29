# Define the grid
grid = [
    [3, 1],
    [3, 2]
]

# Selected positions
selected_positions = [(0, 1), (0, 0), (1, 1)]

# Calculate sums for rows and columns
row_sums = [0, 0]
col_sums = [0, 0]

for r, c in selected_positions:
    row_sums[r] += grid[r][c]
    col_sums[c] += grid[r][c]

# Check if all sums are <= 4
valid = all(s <= 4 for s in row_sums + col_sums)

print(selected_positions if valid else "No valid selection")