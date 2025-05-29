import numpy as np

def count_neighbors(grid, x, y):
    # Count live neighbors with wrapping topology
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

def game_of_life_step(grid):
    rows, cols = grid.shape
    new_grid = np.zeros((rows, cols), dtype=int)
    for x in range(rows):
        for y in range(cols):
            live_neighbors = count_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors in [2, 3]:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
initial_grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,1,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,1,1,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,1,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0,0,0]
])

# Simulate one step
next_grid = game_of_life_step(initial_grid)

# Print the result in JSON format
print(next_grid.tolist())