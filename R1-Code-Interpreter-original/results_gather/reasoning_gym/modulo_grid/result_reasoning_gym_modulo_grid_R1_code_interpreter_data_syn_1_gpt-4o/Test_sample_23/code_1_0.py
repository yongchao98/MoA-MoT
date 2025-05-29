# Initialize a 20x20 grid with all crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with the correct diagonal pattern
for i in range(20):
    grid[i][(i * 4) % 20] = '✅'  # The pattern shifts by 4 columns to the right

# Print the grid
for row in grid:
    print(''.join(row))