# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place ✅ in the 12th column at every second row starting from the second row
for i in range(1, 20, 2):
    grid[i][11] = '✅'

# Convert the grid to a string representation
grid_str = '\n'.join([''.join(row) for row in grid])

print(grid_str)