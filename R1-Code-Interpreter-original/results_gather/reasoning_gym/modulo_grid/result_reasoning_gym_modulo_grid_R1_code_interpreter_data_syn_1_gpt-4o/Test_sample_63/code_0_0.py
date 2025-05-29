# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with ✅ based on the identified pattern
for i in range(17, 20):  # Rows 18 to 20 (0-based index 17 to 19)
    grid[i][i - 17] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))