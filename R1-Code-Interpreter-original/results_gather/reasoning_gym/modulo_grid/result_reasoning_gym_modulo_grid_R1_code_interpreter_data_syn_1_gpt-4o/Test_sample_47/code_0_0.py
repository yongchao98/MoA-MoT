# Initialize a 20x20 grid with "❌"
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the grid with "✅" following the identified pattern
for i in range(20):
    col = (i + 6) % 20
    grid[i][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))