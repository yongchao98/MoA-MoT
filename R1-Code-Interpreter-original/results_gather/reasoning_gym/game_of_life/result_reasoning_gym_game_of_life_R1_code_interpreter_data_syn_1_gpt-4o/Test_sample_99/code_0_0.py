import numpy as np

def get_neighbors(grid, x, y):
    neighbors = [
        ((x - 1) % grid.shape[0], (y - 1) % grid.shape[1]),
        ((x - 1) % grid.shape[0], y),
        ((x - 1) % grid.shape[0], (y + 1) % grid.shape[1]),
        (x, (y - 1) % grid.shape[1]),
        (x, (y + 1) % grid.shape[1]),
        ((x + 1) % grid.shape[0], (y - 1) % grid.shape[1]),
        ((x + 1) % grid.shape[0], y),
        ((x + 1) % grid.shape[0], (y + 1) % grid.shape[1])
    ]
    return neighbors

def game_of_life_step(grid):
    new_grid = np.copy(grid)
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            live_neighbors = sum(grid[nx, ny] for nx, ny in get_neighbors(grid, x, y))
            if grid[x, y] == 1:
                if live_neighbors < 2 or live_neighbors > 3:
                    new_grid[x, y] = 0
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

initial_grid = np.array([
    [0,0,0,0,0,0,1,0,0,0],
    [0,0,1,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1]
])

next_grid = game_of_life_step(initial_grid)
print(next_grid.tolist())