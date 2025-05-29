# Initialize a 20x20 grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with the pattern
for i in range(20):
    grid[i][19 - i] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))