# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(15, 20):  # Rows 16, 17, 18, 19, 20 (0-based index 15, 16, 17, 18, 19)
    grid[i][i - 15] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))