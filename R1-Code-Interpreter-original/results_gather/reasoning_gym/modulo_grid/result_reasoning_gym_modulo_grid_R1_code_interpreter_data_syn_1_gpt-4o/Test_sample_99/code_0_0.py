# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    j = 10 - i
    if 0 <= j < 20:
        grid[i][j] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))