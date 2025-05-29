# Initial grid
grid = [
    [0, 16, 32, 0],
    [0, 32, 8, 0],
    [0, 8, 0, 0],
    [0, 0, 0, 8]
]

# Initial position of the number to move
x, y = 1, 2  # (2, 2) in 0-indexed

# Move sequence
moves = 'urdr'

# Function to perform a move
def move(grid, x, y, direction):
    if direction == 'u' and x > 0:
        if grid[x-1][y] == grid[x][y]:
            grid[x-1][y] *= 2
            grid[x][y] = 0
        elif grid[x-1][y] == 0:
            grid[x-1][y] = grid[x][y]
            grid[x][y] = 0
        x -= 1
    elif direction == 'r' and y < 3:
        if grid[x][y+1] == grid[x][y]:
            grid[x][y+1] *= 2
            grid[x][y] = 0
        elif grid[x][y+1] == 0:
            grid[x][y+1] = grid[x][y]
            grid[x][y] = 0
        y += 1
    elif direction == 'd' and x < 3:
        if grid[x+1][y] == grid[x][y]:
            grid[x+1][y] *= 2
            grid[x][y] = 0
        elif grid[x+1][y] == 0:
            grid[x+1][y] = grid[x][y]
            grid[x][y] = 0
        x += 1
    elif direction == 'l' and y > 0:
        if grid[x][y-1] == grid[x][y]:
            grid[x][y-1] *= 2
            grid[x][y] = 0
        elif grid[x][y-1] == 0:
            grid[x][y-1] = grid[x][y]
            grid[x][y] = 0
        y -= 1
    return x, y

# Apply the moves
for move_direction in moves:
    x, y = move(grid, x, y, move_direction)

# Print the final grid
print(grid)