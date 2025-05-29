# Define the grid
grid = [
    [2, 3],
    [1, 1]
]

# Selected positions
selected_positions = [(0, 1), (1, 0), (1, 1)]

# Calculate the sums for each row and column
row_sums = [0, 0]
col_sums = [0, 0]

for r, c in selected_positions:
    row_sums[r] += grid[r][c]
    col_sums[c] += grid[r][c]

# Check if all sums are less than or equal to 4
valid = all(s <= 4 for s in row_sums + col_sums)

print("Valid selection:", valid)
print("Selected positions:", selected_positions)