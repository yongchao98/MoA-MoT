# Define the grid size
rows, cols = 20, 20

# Initialize the grid with ❌
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Fill in the grid based on the identified pattern
for r in range(rows):
    if r % 4 == 1:  # Rows 2, 6, 10, 14, 18 (0-indexed: 1, 5, 9, 13, 17)
        grid[r][2] = '✅'
        grid[r][17] = '✅'
    elif r % 4 == 3:  # Rows 4, 8, 12, 16, 20 (0-indexed: 3, 7, 11, 15, 19)
        grid[r][8] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))