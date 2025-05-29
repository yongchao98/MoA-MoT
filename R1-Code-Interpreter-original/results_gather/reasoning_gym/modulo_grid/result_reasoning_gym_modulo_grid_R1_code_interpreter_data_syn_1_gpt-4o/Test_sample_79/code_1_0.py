# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid according to the identified pattern
for i in range(20):
    # Place '✅' at the (i+6) % 20 position
    grid[i][(i + 6) % 20] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))