# Initialize the grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for row in range(20):
    col = (7 - row) % 20
    grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))