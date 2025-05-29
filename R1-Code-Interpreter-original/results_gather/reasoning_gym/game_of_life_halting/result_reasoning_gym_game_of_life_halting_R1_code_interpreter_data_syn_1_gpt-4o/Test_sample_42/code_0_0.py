import numpy as np

def get_neighbors_count(grid, x, y):
    # Get the number of live neighbors for a cell at position (x, y)
    neighbors = [
        (-1, -1), (-1, 0), (-1, 1),
        (0, -1),         (0, 1),
        (1, -1), (1, 0), (1, 1)
    ]
    count = 0
    rows, cols = grid.shape
    for dx, dy in neighbors:
        nx, ny = (x + dx) % rows, (y + dy) % cols
        count += grid[nx, ny]
    return count

def game_of_life_step(grid):
    # Perform one step of the Game of Life
    new_grid = np.zeros_like(grid)
    rows, cols = grid.shape
    for x in range(rows):
        for y in range(cols):
            live_neighbors = get_neighbors_count(grid, x, y)
            if grid[x, y] == 1:
                if live_neighbors in [2, 3]:
                    new_grid[x, y] = 1
            else:
                if live_neighbors == 3:
                    new_grid[x, y] = 1
    return new_grid

def will_halt_before_20_steps(initial_grid):
    grid = np.array(initial_grid)
    for step in range(20):
        grid = game_of_life_step(grid)
        if np.all(grid == 0):
            return True
    return False

initial_board = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0],
    [0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

result = will_halt_before_20_steps(initial_board)
print(result)