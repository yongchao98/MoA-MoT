# Define the grid size
grid_size = 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Define the positions of '✅' based on the observed pattern
positions = [
    (1, 8), (1, 18),
    (3, 2), (3, 13),
    (4, 1),  # This is the position of the '❔' which should be '✅'
    (7, 6), (7, 17),
    (9, 7), (9, 18),
    (11, 8), (11, 19),
    (13, 2), (13, 13),
    (17, 6), (17, 17),
    (19, 7), (19, 18)
]

# Place '✅' in the grid at the specified positions
for x, y in positions:
    grid[x][y] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))