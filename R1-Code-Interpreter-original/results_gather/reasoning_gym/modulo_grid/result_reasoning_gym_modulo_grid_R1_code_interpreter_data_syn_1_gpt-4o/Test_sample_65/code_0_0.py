# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with the pattern
for i in range(20):
    if i < 10:
        grid[i][9 - i] = '✅'
    else:
        grid[i][i - 10] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))