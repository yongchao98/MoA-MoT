# Define the grid size
rows, cols = 20, 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Define the pattern for '✅'
pattern_rows_1 = [3, 9, 15]
pattern_cols_1 = [9, 17]

pattern_rows_2 = [5, 11, 17]
pattern_cols_2 = [3, 6, 9, 12, 15, 18]

# Fill the grid according to the pattern
for r in pattern_rows_1:
    for c in pattern_cols_1:
        grid[r][c] = '✅'

for r in pattern_rows_2:
    for c in pattern_cols_2:
        grid[r][c] = '✅'

# Convert the grid to a string format
completed_grid = '\n'.join([''.join(row) for row in grid])

# Print the completed grid
print(completed_grid)