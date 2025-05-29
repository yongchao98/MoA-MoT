# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Define the pattern for the columns where ✅ appears
pattern_columns = [2, 8, 12, 18]  # 0-indexed

# Fill in the grid according to the pattern
for row in range(1, 20, 2):  # 1-indexed rows: 2, 4, 6, ..., 20
    for col in pattern_columns:
        grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))