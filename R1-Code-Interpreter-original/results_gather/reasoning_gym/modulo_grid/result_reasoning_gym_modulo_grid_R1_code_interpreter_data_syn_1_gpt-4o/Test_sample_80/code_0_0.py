# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified patterns
for row in range(20):
    for col in range(20):
        if row + col == 14 or row - col == 14:
            grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))