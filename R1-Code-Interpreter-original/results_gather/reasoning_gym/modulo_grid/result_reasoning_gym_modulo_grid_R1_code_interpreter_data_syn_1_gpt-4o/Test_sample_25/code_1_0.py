# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place the known check mark
grid[1][10] = '✅'

# Fill the grid based on the vertical pattern hypothesis
for i in range(20):
    grid[(1 + i) % 20][10] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))