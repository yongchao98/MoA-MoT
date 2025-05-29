# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the corrected pattern
for i in range(20):
    if i < 10:
        grid[i][i + 10] = '✅'
    else:
        grid[i][i - 10] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))