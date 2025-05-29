import numpy as np

def count_neighbors(grid, x, y):
    # Count live neighbors with wrapping
    neighbors = [
        ((x-1) % grid.shape[0], (y-1) % grid.shape[1]),
        ((x-1) % grid.shape[0], y),
        ((x-1) % grid.shape[0], (y+1) % grid.shape[1]),
        (x, (y-1) % grid.shape[1]),
        (x, (y+1) % grid.shape[1]),
        ((x+1) % grid.shape[0], (y-1) % grid.shape[1]),
        ((x+1) % grid.shape[0], y),
        ((x+1) % grid.shape[0], (y+1) % grid.shape[1])
    ]
    return sum(grid[nx, ny] for nx, ny in neighbors)

def game_of_life_step(grid):
    new_grid = np.zeros(grid.shape, dtype=int)
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            live_neighbors = count_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors in [2, 3]:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,1,0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,1,0,0,0],
    [0,0,0,0,0,0,0,1,1,0],
    [0,1,1,0,0,0,0,0,0,0]
])

# Simulate one step
new_grid = game_of_life_step(grid)

# Print the new grid in JSON format
print(new_grid.tolist())