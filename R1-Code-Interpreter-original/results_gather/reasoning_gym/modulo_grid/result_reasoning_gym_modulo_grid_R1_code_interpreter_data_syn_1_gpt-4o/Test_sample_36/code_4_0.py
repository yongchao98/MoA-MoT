# Define the grid size
grid_size = 20

# Initialize the grid with False
grid = [[False for _ in range(grid_size)] for _ in range(grid_size)]

# Mark the known checkmarks
grid[19][0] = True
grid[0][19] = True

# Simulate a path from (19, 0) to (0, 19)
x, y = 19, 0
while x >= 0 and y < grid_size:
    grid[x][y] = True
    x -= 1
    y += 1

# Check if the position (13, 15) is marked
is_marked = grid[13][15]

print(is_marked)