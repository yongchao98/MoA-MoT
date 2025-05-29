# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the grid based on the identified pattern
for i in range(20):
    if i <= 8:
        grid[i][8 + i] = '✅'
    if i >= 8:
        grid[i][i - 8] = '✅'
        if i + 8 < 20:
            grid[i][i + 8] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))