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
    new_grid = np.copy(grid)
    # Iterate over each cell in the grid
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            # Get the neighbors
            neighbors = get_neighbors(grid, x, y)
            # Count the number of live neighbors
            live_neighbors = sum(grid[nx, ny] for nx, ny in neighbors)
            # Apply the Game of Life rules
            if grid[x, y] == 1:
                if live_neighbors < 2 or live_neighbors > 3:
                    new_grid[x, y] = 0
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,1,0,0,0],
    [0,0,1,0,0,0,0,0,0,1],
    [0,0,0,0,1,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [1,0,0,1,0,0,1,0,0,0]
])

# Simulate one step
new_grid = game_of_life_step(grid)

# Convert the result to a list of lists and print it
print(new_grid.tolist())