# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    for j in range(20):
        if i + j == 4 + i:
            grid[i][j] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))