# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the first half of the grid
for i in range(10):
    grid[i][19-i] = '✅'

# Fill the second half of the grid
for i in range(10, 20):
    grid[i][i-10] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))