# Initialize a 20x20 grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the grid based on the identified pattern
for i in range(20):
    grid[i][19 - i] = '✅'  # Top right to bottom left diagonal
    grid[i][i] = '✅'       # Bottom left to top right diagonal

# Print the grid
for row in grid:
    print(''.join(row))