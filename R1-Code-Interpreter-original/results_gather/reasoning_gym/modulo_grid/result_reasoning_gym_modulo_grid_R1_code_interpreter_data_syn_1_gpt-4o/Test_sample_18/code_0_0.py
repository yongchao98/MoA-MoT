# Define the grid size
grid_size = 20

# Initialize the grid with crosses
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Define the positions of the checkmarks
checkmark_positions = [
    (0, 9), (3, 3), (3, 19), (5, 5), (7, 15), (9, 1), (9, 17),
    (11, 11), (13, 13), (15, 7), (17, 9), (19, 3), (19, 19)
]

# Place the checkmarks in the grid
for row, col in checkmark_positions:
    grid[row][col] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))