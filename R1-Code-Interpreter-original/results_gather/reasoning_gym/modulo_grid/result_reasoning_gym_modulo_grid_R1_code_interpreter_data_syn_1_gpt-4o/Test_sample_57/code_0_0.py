# Define the grid size
rows, cols = 20, 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Place '✅' according to the identified pattern
for row in range(4, rows, 5):
    grid[row][5] = '✅'
    grid[row][12] = '✅'

# Replace the '?' with '❌'
grid[1][16] = '❌'

# Print the completed grid
for row in grid:
    print(''.join(row))