# Define the grid size
rows, cols = 20, 20

# Initialize the grid with ❌
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Fill in the pattern with ✅
for row in range(1, rows, 2):
    for col in [2, 8, 14]:
        grid[row][col] = '✅'

# Convert the grid to a string representation
completed_grid = '\n'.join([''.join(row) for row in grid])

# Print the completed grid
print(completed_grid)