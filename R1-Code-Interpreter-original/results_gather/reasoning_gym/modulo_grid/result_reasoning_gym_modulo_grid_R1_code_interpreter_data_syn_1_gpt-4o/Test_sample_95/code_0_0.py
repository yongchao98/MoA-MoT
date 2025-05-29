# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for row in range(20):
    for col in range(20):
        if row % 5 == col % 5:
            grid[row][col] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))