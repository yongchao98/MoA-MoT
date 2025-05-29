# Initialize the grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the check marks according to the identified pattern
for i in range(20):
    grid[i][(i + 10) % 20] = '✅'

# Convert the grid to a string representation
completed_grid = '\n'.join([''.join(row) for row in grid])

print(completed_grid)