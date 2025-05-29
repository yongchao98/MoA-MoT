# Initialize the grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with check marks based on the identified pattern
for row in range(20):
    col = (row + 2) % 19
    grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))