# Define the grid size
rows, cols = 20, 20

# Initialize the grid with crosses
grid = [['❌' for _ in range(cols)] for _ in range(rows)]

# Fill in the check marks based on the identified pattern
for r in range(rows):
    if r % 4 == 1:  # Rows 2, 6, 10, 14, 18 (0-indexed: 1, 5, 9, 13, 17)
        for c in range(3, cols, 4):  # Columns 4, 8, 12, 16, 20 (0-indexed: 3, 7, 11, 15, 19)
            grid[r][c] = '✅'
    elif r % 4 == 3:  # Rows 4, 8, 12, 16, 20 (0-indexed: 3, 7, 11, 15, 19)
        for c in range(1, cols, 4):  # Columns 2, 6, 10, 14, 18 (0-indexed: 1, 5, 9, 13, 17)
            grid[r][c] = '✅'

# Convert the grid to a string representation
completed_grid = '\n'.join(''.join(row) for row in grid)
print(completed_grid)