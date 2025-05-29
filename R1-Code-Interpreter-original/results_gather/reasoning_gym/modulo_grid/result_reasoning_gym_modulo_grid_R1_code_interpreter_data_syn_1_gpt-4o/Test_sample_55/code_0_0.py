# Define the grid size
grid_size = 20

# Initialize the grid with '❌'
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Define the positions of '✅' based on the identified pattern
positions = [
    (1, 3),
    (3, 1), (3, 7), (3, 13), (3, 19),
    (5, 15),
    (7, 3),
    (11, 15),
    (13, 3),
    (15, 5), (15, 11), (15, 17),
    (17, 15),
    (19, 3)
]

# Place '✅' in the grid
for x, y in positions:
    grid[x][y] = '✅'

# Convert the grid to a string for display
grid_str = '\n'.join(''.join(row) for row in grid)
print(grid_str)