# Define the grid size
rows, cols = 20, 20

# Initialize the grid with ❌
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Apply the vertical pattern
for r in range(1, rows, 5):
    for c in range(1, cols, 4):
        grid[r][c] = '✅'

# Apply the horizontal pattern
for r in range(3, rows, 4):
    for c in range(3, cols, 4):
        grid[r][c] = '✅'

# Convert the grid to a string format
completed_grid = '\n'.join([''.join(row) for row in grid])

print(completed_grid)