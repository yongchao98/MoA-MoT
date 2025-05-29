import numpy as np

def get_live_neighbors_count(grid, x, y):
    # Get the size of the grid
    rows, cols = grid.shape
    # Calculate the sum of the 3x3 block centered around (x, y)
    # Use modulo to wrap around the edges
    live_neighbors = sum(grid[(x + i) % rows, (y + j) % cols] 
                         for i in [-1, 0, 1] for j in [-1, 0, 1]) - grid[x, y]
    return live_neighbors

def game_of_life_step(grid):
    # Create a copy of the grid to store the next state
    next_grid = np.copy(grid)
    # Iterate over each cell in the grid
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            live_neighbors = get_live_neighbors_count(grid, x, y)
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
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,1,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,1,0,0,1,0],
    [0,1,0,0,0,0,0,0,0,0]
])

# Simulate one step
next_grid = game_of_life_step(initial_grid)

# Convert the result to a list of lists and print it
print(next_grid.tolist())