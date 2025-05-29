import numpy as np

# Initial grid
grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,1,1,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,1,0,0,1,0,0,0,0],
    [0,1,0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
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
    rows, cols = grid.shape
    for dx, dy in neighbors:
        nx, ny = (x + dx) % rows, (y + dy) % cols
        count += grid[nx, ny]
    return count

# Function to simulate one step
def simulate_step(grid):
    rows, cols = grid.shape
    new_grid = np.zeros((rows, cols), dtype=int)
    for x in range(rows):
        for y in range(cols):
            live_neighbors = count_live_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors == 2 or live_neighbors == 3:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Simulate one step
new_grid = simulate_step(grid)

# Convert to list of lists for JSON format
result = new_grid.tolist()
print(result)