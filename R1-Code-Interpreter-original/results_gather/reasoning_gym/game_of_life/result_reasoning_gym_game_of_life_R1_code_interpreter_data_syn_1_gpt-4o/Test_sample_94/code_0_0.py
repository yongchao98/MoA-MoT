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

def game_of_life_step(grid):
    """Simulate one step of the Game of Life."""
    new_grid = np.zeros(grid.shape, dtype=int)
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            live_neighbors = get_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors in (2, 3):
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,1,0,1,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,1,0,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,1]
])

# Simulate one step
new_grid = game_of_life_step(grid)

# Convert to list of lists for JSON format
result = new_grid.tolist()
print(result)