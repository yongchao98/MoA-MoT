# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Place '✅' according to the identified pattern
for i in range(20):
    grid[i][(i + 6) % 20] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))