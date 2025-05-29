import numpy as np

# Initial grid
grid = np.array([
    [1,0,0,0,0,0,0,0,1,0],
    [0,0,0,1,0,1,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,1,0,0,1,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,0]
])

# Function to count live neighbors with wrapping
def count_live_neighbors(grid, x, y):
    neighbors = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),         (0, 1),
        (1, -1), (1, 0), (1, 1)
    ]
    count = 0
    for dx, dy in neighbors:
        nx, ny = (x + dx) % grid.shape[0], (y + dy) % grid.shape[1]
        count += grid[nx, ny]
    return count

# Create a new grid for the next state
new_grid = np.zeros_like(grid)

# Apply the rules of the Game of Life
for x in range(grid.shape[0]):
    for y in range(grid.shape[1]):
        live_neighbors = count_live_neighbors(grid, x, y)
        if grid[x, y] == 1:
            if live_neighbors in [2, 3]:
                new_grid[x, y] = 1
        else:
            if live_neighbors == 3:
                new_grid[x, y] = 1

# Convert the result to a list of lists and print
result = new_grid.tolist()
print(result)