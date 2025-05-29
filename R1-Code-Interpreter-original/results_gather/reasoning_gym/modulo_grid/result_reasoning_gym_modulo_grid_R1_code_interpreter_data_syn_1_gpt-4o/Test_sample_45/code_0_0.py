# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' based on the identified pattern
for i in range(20):
    col = 8 - (i % 8) - 1  # Calculate the column index for '✅'
    grid[i][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))