import numpy as np

def get_neighbors(grid, x, y):
    neighbors = []
    rows, cols = grid.shape
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx != 0 or dy != 0:
                nx, ny = (x + dx) % rows, (y + dy) % cols
                neighbors.append(grid[nx, ny])
    return neighbors

def game_of_life_step(grid):
    rows, cols = grid.shape
    new_grid = np.zeros((rows, cols), dtype=int)
    for x in range(rows):
        for y in range(cols):
            live_neighbors = sum(get_neighbors(grid, x, y))
            if grid[x, y] == 1:
                if live_neighbors in [2, 3]:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

initial_grid = np.array([
    [0,0,0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,1,1,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [1,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,1,0,0,0,0,0,0,0]
])

next_grid = game_of_life_step(initial_grid)
print(next_grid.tolist())