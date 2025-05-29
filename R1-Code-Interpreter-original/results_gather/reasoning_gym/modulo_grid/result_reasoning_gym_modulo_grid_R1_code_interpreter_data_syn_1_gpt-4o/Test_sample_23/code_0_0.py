# Initialize a 20x20 grid with all crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with the diagonal pattern
for i in range(20):
    grid[i][(i + 3) % 20] = '✅'  # The pattern shifts by 3 columns to the right

# Print the grid
for row in grid:
    print(''.join(row))