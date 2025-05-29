import numpy as np

def count_neighbors(grid, x, y):
    # Get the size of the grid
    rows, cols = grid.shape
    # Calculate the sum of the 3x3 block centered on (x, y)
    # Use modulo to wrap around the edges
    return sum(grid[(x + i) % rows, (y + j) % cols] for i in range(-1, 2) for j in range(-1, 2)) - grid[x, y]

def game_of_life_step(grid):
    # Create a copy of the grid to store the next state
    new_grid = grid.copy()
    # Iterate over each cell in the grid
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            # Count the number of live neighbors
            live_neighbors = count_neighbors(grid, x, y)
            # Apply the rules of the Game of Life
            if grid[x, y] == 1:
                if live_neighbors < 2 or live_neighbors > 3:
                    new_grid[x, y] = 0
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
initial_grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,1,0],
    [1,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,1,0,1,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0]
])

# Simulate one step
next_grid = game_of_life_step(initial_grid)

# Convert the result to a list of lists and print it
print(next_grid.tolist())