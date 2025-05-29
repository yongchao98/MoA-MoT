# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' marks according to the identified pattern
for row in range(20):
    col = (row + 10) % 20  # Calculate the column index for '✅'
    grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))