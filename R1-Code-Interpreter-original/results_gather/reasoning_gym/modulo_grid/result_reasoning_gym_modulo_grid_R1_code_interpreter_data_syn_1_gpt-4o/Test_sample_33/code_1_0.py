# Define the grid size
grid_size = 20

# Define the positions of the "✅" marks, including the missing one
positions = [
    (2, 11), (3, 6), (3, 13), (3, 20), (4, 9), (5, 7), (5, 14),
    (6, 3), (6, 17), (7, 5), (7, 12), (7, 19), (9, 4), (9, 10),
    (9, 17), (10, 13), (11, 2), (11, 9), (11, 16), (12, 7),
    (13, 3), (13, 10), (13, 17), (14, 5), (14, 19), (16, 11),
    (17, 6), (17, 8), (17, 13), (17, 20), (18, 9), (19, 7),
    (19, 14), (20, 3), (20, 17)
]

# Initialize the grid with "❌"
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Fill the grid with "✅" at the specified positions
for row, col in positions:
    grid[row - 1][col - 1] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))