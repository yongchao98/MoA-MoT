# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    if i + 6 < 20:
        grid[i][i + 6] = '✅'
    if i - 14 >= 0:
        grid[i][i - 14] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))