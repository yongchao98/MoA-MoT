# Initialize a 20x20 grid with all ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the grid based on the identified pattern
for i in range(20):
    for j in range(20):
        if i + j == 6:  # The sum of indices for the checkmark pattern
            grid[i][j] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))