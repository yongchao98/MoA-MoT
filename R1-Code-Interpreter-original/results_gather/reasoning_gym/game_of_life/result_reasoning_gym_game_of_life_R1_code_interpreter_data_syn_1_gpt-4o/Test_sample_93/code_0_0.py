import numpy as np

def count_neighbors(grid, x, y):
    # Count live neighbors in a wrapping grid
    neighbors = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),         (0, 1),
        (1, -1), (1, 0), (1, 1)
    ]
    count = 0
    for dx, dy in neighbors:
        nx, ny = (x + dx) % grid.shape[0], (y + dy) % grid.shape[1]
        count += grid[nx, ny]
    return count

def game_of_life_step(grid):
    new_grid = np.zeros_like(grid)
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            live_neighbors = count_neighbors(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors == 2 or live_neighbors == 3:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

# Initial grid
grid = np.array([
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,0,0],
    [1,0,1,0,1,0,0,0,0,0]
])

# Simulate one step
new_grid = game_of_life_step(grid)

# Convert to list of lists for JSON format
result = new_grid.tolist()
print(result)