# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with '✅' based on the identified pattern
for i in range(20):
    if i < 13:
        grid[i][12 - i] = '✅'
    else:
        grid[i][i - 7] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))