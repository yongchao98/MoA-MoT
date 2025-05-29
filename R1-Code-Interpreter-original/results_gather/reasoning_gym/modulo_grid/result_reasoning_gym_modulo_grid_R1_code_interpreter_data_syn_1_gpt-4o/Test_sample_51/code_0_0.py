# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with ✅ based on the identified pattern
for i in range(20):
    j = 15 - i
    if 0 <= j < 20:
        grid[i][j] = '✅'

# Convert the grid to a string format for display
completed_grid = '\n'.join([''.join(row) for row in grid])
print(completed_grid)