# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place ✅ at positions (i, i) for i from 0 to 19
for i in range(20):
    grid[i][i] = '✅'

# Convert the grid to a string format for display
completed_grid = '\n'.join([''.join(row) for row in grid])

print(completed_grid)