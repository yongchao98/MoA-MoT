# Initialize the grid with crosses
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    j = (4 - i) % 20  # Calculate the column index for the checkmark
    grid[i][j] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))