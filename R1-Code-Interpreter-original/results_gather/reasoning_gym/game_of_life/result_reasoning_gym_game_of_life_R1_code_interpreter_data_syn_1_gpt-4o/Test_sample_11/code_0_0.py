import numpy as np

def count_neighbors(grid, x, y):
    # Count live neighbors with wrapping
    rows, cols = grid.shape
    return sum(grid[(x + i) % rows, (y + j) % cols] 
               for i in [-1, 0, 1] 
               for j in [-1, 0, 1] 
               if not (i == 0 and j == 0))

def game_of_life_step(grid):
    rows, cols = grid.shape
    new_grid = np.zeros((rows, cols), dtype=int)
    
    for x in range(rows):
        for y in range(cols):
            live_neighbors = count_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors in [2, 3]:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
initial_grid = np.array([
    [0,0,0,0,1,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,1],
    [0,0,1,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,0,0,0,1],
    [0,0,0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,0,0,0,0]
])

# Simulate one step
next_grid = game_of_life_step(initial_grid)

# Convert to list of lists for JSON format
next_grid_list = next_grid.tolist()

print(next_grid_list)