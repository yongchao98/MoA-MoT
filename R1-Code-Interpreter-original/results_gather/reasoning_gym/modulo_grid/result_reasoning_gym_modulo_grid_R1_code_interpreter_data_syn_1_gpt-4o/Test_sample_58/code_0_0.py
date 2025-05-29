# Initialize a 20x20 grid with "❌"
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the "✅" according to the identified pattern
for row in range(0, 20, 6):
    for col in range(3, 20, 6):
        grid[row][col] = '✅'

for row in range(5, 20, 6):
    for col in range(5, 20, 6):
        grid[row][col] = '✅'

# Convert the grid to a string format
completed_grid = '\n'.join([''.join(row) for row in grid])
print(completed_grid)