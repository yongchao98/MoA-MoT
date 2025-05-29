# Define the grid size
grid_size = 20

# Initialize the grid with ❌
grid = [['❌' for _ in range(grid_size)] for _ in range(grid_size)]

# Define the pattern for the first 10 rows
pattern = [
    (0, 10),
    (1, 6), (1, 7), (1, 19),
    (4, 4), (4, 17),
    (6, 10),
    (9, 2), (9, 11), (9, 15)
]

# Apply the pattern to the grid
for row in range(0, grid_size, 10):
    for (r, c) in pattern:
        grid[(r + row) % grid_size][c] = '✅'

# Fill the missing symbol at (18, 18)
grid[18][18] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))