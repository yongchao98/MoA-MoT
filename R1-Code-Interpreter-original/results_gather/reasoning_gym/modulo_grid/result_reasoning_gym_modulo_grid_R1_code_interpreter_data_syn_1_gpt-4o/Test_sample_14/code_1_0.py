# Initialize the grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the check marks according to the identified pattern
for i in range(6):
    grid[i][5 - i] = '✅'
for i in range(6, 20):
    grid[i][23 - i] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))