import numpy as np

def get_wrapped_index(index, size):
    return index % size

def count_live_neighbors(grid, x, y):
    directions = [(-1, -1), (-1, 0), (-1, 1), 
                  (0, -1),         (0, 1), 
                  (1, -1), (1, 0), (1, 1)]
    live_neighbors = 0
    rows, cols = grid.shape
    for dx, dy in directions:
        nx, ny = get_wrapped_index(x + dx, rows), get_wrapped_index(y + dy, cols)
        live_neighbors += grid[nx, ny]
    return live_neighbors

def game_of_life_step(grid):
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

initial_grid = np.array([[0,0,0,0,0,0,0,0,0,1],
                         [0,0,0,0,0,1,0,0,0,0],
                         [0,0,0,0,0,0,0,0,0,0],
                         [1,0,0,0,0,0,0,0,0,0],
                         [0,0,0,0,0,0,0,0,1,0],
                         [0,0,0,0,1,0,0,0,0,0],
                         [0,0,0,1,0,0,0,0,0,0],
                         [0,0,0,0,0,0,0,0,1,0],
                         [1,0,0,0,0,0,0,0,0,0],
                         [0,1,0,1,0,0,0,0,0,0]])

next_grid = game_of_life_step(initial_grid)
print(next_grid.tolist())