# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid with '✅' based on the identified pattern
for i in range(20):
    grid[i][19 - i] = '✅'  # Diagonal from top-left to bottom-right
    grid[i][i] = '✅'  # Diagonal from bottom-left to top-right

# Print the completed grid
for row in grid:
    print(''.join(row))