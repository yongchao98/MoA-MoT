# Define the grid size
grid_size = 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Fill the grid with '✅' following the identified pattern
for i in range(grid_size):
    position = (i + 7) % grid_size  # Calculate the position of '✅' in each row
    grid[i][position] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))