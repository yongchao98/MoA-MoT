# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with '✅' based on the identified pattern
for i in range(20):
    grid[i][(i + 9) % 20] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))