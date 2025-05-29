# Define the grid size
grid_size = 20

# Initialize the grid with crosses
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Place the checkmarks
for row, col in checkmarks:
    grid[row - 1][col - 1] = '✅'

# Place the checkmark for the question mark position
grid[3 - 1][11 - 1] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))