# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with '✅' according to the identified pattern
for i in range(20):
    grid[i][19 - i] = '✅'  # Top right to bottom left diagonal
    grid[i][i] = '✅'       # Top left to bottom right diagonal

# Convert the grid to a string format for display
completed_grid = '\n'.join([''.join(row) for row in grid])
print(completed_grid)