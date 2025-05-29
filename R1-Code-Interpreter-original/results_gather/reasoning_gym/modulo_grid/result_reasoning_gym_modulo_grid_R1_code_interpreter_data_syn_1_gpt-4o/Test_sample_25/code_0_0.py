# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place the known check mark
grid[1][10] = '✅'

# Fill the grid based on the diagonal pattern hypothesis
for i in range(20):
    if 0 <= i + 9 < 20:  # Ensure we stay within bounds
        grid[i][i + 9] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))