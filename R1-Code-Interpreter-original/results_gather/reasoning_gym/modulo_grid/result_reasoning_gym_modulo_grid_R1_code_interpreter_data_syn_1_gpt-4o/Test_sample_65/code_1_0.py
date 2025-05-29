# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with the pattern
for i in range(20):
    if i < 10:
        grid[i][i + 2] = '✅'
    else:
        grid[i][19 - (i - 10)] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))