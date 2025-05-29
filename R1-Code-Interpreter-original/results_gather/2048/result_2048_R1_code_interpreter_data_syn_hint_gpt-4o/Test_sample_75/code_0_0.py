# Initial grid
grid = [
    [0, 16, 4, 16, 0],
    [0, 2, 0, 0, 32],
    [0, 0, 0, 0, 16],
    [4, 0, 0, 0, 8],
    [0, 0, 0, 2, 0]
]

# Initial position
x, y = 3, 3

# Move sequence
moves = 'urrld'

# Move directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Simulate the moves
for move in moves:
    dx, dy = directions[move]
    x_new, y_new = x + dx, y + dy
    # Check if the move is within bounds
    if 0 <= x_new < 5 and 0 <= y_new < 5:
        # Check if the numbers can combine
        if grid[x_new][y_new] == grid[x][y] and grid[x][y] != 0:
            grid[x_new][y_new] *= 2
            grid[x][y] = 0
        elif grid[x][y] != 0:
            grid[x_new][y_new] = grid[x][y]
            grid[x][y] = 0
        x, y = x_new, y_new

# Print the final grid
print(grid)