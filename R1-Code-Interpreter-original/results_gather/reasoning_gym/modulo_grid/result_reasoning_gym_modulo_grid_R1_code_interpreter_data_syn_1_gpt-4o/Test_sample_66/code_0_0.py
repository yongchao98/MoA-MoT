# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' symbols based on the identified pattern
for i in range(20):
    grid[19 - i][i] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))