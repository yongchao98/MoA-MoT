# Define the grid size
grid_size = 20

# Initialize the grid with ❌
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Fill in the check marks (✅) following the identified pattern
for i in range(grid_size):
    grid[i][grid_size - i - 1] = '✅'

# Fill in the missing check mark (❔) in the 5th row, 7th column
grid[4][6] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))