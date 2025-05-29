# Initialize the grid
grid_size = 20
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Fill the grid based on the identified pattern
for i in range(grid_size):
    grid[i][i % grid_size] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))