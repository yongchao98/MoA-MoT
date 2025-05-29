# Define the grid with the pattern
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid based on the identified pattern
for i in range(20):
    if i % 4 == 1 or i % 4 == 3:  # Rows 2, 6, 8, 12, 14, 18, 20
        grid[i][7] = '✅'
        grid[i][16] = '✅'
    elif i % 4 == 2:  # Rows 4, 10, 16
        for j in range(2, 19, 3):
            grid[i][j] = '✅'

# Convert the grid to a string representation
grid_str = '\n'.join([''.join(row) for row in grid])
print(grid_str)