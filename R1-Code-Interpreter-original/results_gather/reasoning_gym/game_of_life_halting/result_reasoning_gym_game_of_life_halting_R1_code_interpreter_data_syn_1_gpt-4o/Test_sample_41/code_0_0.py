import numpy as np

def get_neighbors(grid, x, y):
    """Get the number of live neighbors for a cell at position (x, y) with wrapping."""
    neighbors = 0
    rows, cols = grid.shape
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                continue
            neighbors += grid[(x + i) % rows, (y + j) % cols]
    return neighbors

def step(grid):
    """Perform one step in the Game of Life."""
    new_grid = np.zeros_like(grid)
    rows, cols = grid.shape
    for x in range(rows):
        for y in range(cols):
            live_neighbors = get_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors == 2 or live_neighbors == 3:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

def will_halt(initial_grid, max_steps=20):
    """Check if the game will halt at or before max_steps."""
    grid = np.array(initial_grid)
    for _ in range(max_steps):
        grid = step(grid)
        if np.all(grid == 0):
            return True
    return False

initial_board = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

print(will_halt(initial_board))