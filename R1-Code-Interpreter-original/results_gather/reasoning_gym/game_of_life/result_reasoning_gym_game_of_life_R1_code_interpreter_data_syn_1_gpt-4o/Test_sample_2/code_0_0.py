import numpy as np

def get_neighbors(grid, x, y):
    # Get the size of the grid
    rows, cols = grid.shape
    # Calculate the indices of the neighbors with wrapping
    neighbors = [
        ((x - 1) % rows, (y - 1) % cols), ((x - 1) % rows, y), ((x - 1) % rows, (y + 1) % cols),
        (x, (y - 1) % cols),                             (x, (y + 1) % cols),
        ((x + 1) % rows, (y - 1) % cols), ((x + 1) % rows, y), ((x + 1) % rows, (y + 1) % cols)
    ]
    return neighbors

def game_of_life_step(grid):
    # Create a copy of the grid to store the next state
    next_grid = np.copy(grid)
    # Iterate over each cell in the grid
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            # Count the number of live neighbors
            live_neighbors = sum(grid[nx, ny] for nx, ny in get_neighbors(grid, x, y))
            # Apply the rules of the Game of Life
            if grid[x, y] == 1:
                if live_neighbors < 2 or live_neighbors > 3:
                    next_grid[x, y] = 0
            else:
                if live_neighbors == 3:
                    next_grid[x, y] = 1
    return next_grid

# Initial grid
initial_grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,1,0,0,1],
    [0,0,0,0,0,0,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,1,0,0,0,0],
    [1,0,0,0,0,0,0,0,0,0]
])

# Simulate one step
next_grid = game_of_life_step(initial_grid)

# Convert the result to a list of lists and print it
result = next_grid.tolist()
print(result)