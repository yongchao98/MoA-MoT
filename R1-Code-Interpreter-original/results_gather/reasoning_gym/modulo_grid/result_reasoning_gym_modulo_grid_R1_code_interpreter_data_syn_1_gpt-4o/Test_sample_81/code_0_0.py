# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place '✅' at every 9th row in the 10th column (0-based index)
for i in range(1, 20, 9):
    grid[i][10] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))