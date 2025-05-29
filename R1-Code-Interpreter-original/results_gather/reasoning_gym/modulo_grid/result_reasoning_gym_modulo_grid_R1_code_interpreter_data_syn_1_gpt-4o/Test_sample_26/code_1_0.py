# Initialize the grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the checkmarks based on the identified pattern
for i in range(20):
    if i + 13 < 20:
        grid[i][i + 13] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))