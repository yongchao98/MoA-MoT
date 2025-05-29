import numpy as np

def get_wrapped_index(index, max_index):
    """Get the wrapped index for a toroidal array."""
    return index % max_index

def count_live_neighbors(grid, x, y):
    """Count the number of live neighbors for a cell at position (x, y)."""
    rows, cols = grid.shape
    live_neighbors = 0
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx != 0 or dy != 0:  # Exclude the cell itself
                nx = get_wrapped_index(x + dx, rows)
                ny = get_wrapped_index(y + dy, cols)
                live_neighbors += grid[nx, ny]
    return live_neighbors

def game_of_life_step(grid):
    """Perform one step of the Game of Life."""
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

# Initial grid
initial_grid = np.array([
    [0,0,0,0,0,1,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,0,1],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,1,1,0,1,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0]
])

# Perform one step
new_grid = game_of_life_step(initial_grid)

# Convert to list of lists for JSON format
new_grid_list = new_grid.tolist()

print(new_grid_list)