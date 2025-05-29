# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern for '✅'
pattern_rows = [2, 8, 14, 20]  # 0-indexed rows for pattern 1
pattern_cols_1 = [8, 16]       # 0-indexed columns for pattern 1

pattern_rows_2 = [4, 10, 16]   # 0-indexed rows for pattern 2
pattern_cols_2 = [2, 5, 8, 11, 14, 17]  # 0-indexed columns for pattern 2

# Apply the pattern to the grid
for row in pattern_rows:
    for col in pattern_cols_1:
        grid[row][col] = '✅'

for row in pattern_rows_2:
    for col in pattern_cols_2:
        grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))