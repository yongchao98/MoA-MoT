# Define the grid size
grid_size = 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# List of positions for '✅'
positions = [
    (1, 18), (2, 9), (3, 4), (3, 6), (4, 14), (5, 15), (6, 16),
    (7, 3), (8, 8), (9, 7), (10, 2), (11, 17), (12, 12), (13, 11),
    (14, 16), (15, 5), (16, 6), (17, 13), (18, 10), (19, 1)
]

# Fill the grid with '✅' at the specified positions
for row, col in positions:
    grid[row - 1][col - 1] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))