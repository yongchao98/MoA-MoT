# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with the pattern
for i in range(20):
    grid[i][19 - i] = '✅'

# Replace the question mark with the correct pattern
grid[7][14] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))