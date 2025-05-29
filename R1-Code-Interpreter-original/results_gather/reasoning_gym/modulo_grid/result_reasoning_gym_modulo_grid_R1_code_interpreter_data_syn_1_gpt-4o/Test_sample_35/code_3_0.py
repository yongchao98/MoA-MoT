# Define the grid size
grid_size = 20

# Initialize the grid with "❌"
grid = [["❌" for _ in range(grid_size)] for _ in range(grid_size)]

# Define the positions of the "✅" symbols
positions = [
    (1, 14), (2, 8), (4, 11), (8, 3), (8, 17), 
    (9, 13), (12, 5), (12, 19), (13, 9), (14, 8), 
    (15, 2), (15, 16), (17, 14), (18, 8), (20, 12)
]

# Place the "✅" symbols in the grid
for row, col in positions:
    grid[row-1][col-1] = "✅"

# Print the grid
for row in grid:
    print("".join(row))