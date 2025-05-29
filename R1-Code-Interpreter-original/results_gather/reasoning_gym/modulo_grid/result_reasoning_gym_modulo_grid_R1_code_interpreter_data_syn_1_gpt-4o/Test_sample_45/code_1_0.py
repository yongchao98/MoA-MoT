# Initialize the grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' based on the corrected pattern
for i in range(20):
    grid[i][(7 - i % 8)] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))