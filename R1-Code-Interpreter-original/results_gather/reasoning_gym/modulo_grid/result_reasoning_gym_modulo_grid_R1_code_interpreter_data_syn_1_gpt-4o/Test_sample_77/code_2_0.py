# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern for '✅'
pattern_cols_1 = [8, 16]     # 0-indexed columns for pattern 1
pattern_cols_2 = [2, 5, 8, 11, 14, 17]  # 0-indexed columns for pattern 2

# Apply the pattern to the grid
for row in range(20):
    if row % 6 == 2:  # Pattern 1 rows
        for col in pattern_cols_1:
            grid[row][col] = '✅'
    elif row % 6 == 3:  # Pattern 2 rows
        for col in pattern_cols_2:
            grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))