# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the checkmarks based on the identified pattern
for row in range(1, 20, 2):  # Every even row (1-indexed)
    grid[row][3] = '✅'  # 4th column (0-indexed)
    grid[row][10] = '✅'  # 11th column (0-indexed)
    grid[row][19] = '✅'  # 20th column (0-indexed)

# Convert the grid to a string representation
completed_grid = '\n'.join([''.join(row) for row in grid])
print(completed_grid)