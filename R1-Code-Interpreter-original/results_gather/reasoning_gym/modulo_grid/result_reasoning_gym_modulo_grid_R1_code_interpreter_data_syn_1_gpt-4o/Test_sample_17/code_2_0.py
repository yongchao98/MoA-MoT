# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with '✅' according to the identified pattern
for i in range(20):
    col = (i + 2) % 20  # Adjusted pattern to match the original grid
    grid[i][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))