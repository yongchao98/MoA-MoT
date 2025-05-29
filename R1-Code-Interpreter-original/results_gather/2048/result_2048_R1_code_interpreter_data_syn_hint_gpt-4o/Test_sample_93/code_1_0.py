# Initial grid
grid = [
    [0, 0, 0, 0, 2],
    [0, 32, 2, 0, 0],
    [0, 4, 16, 0, 2],
    [32, 0, 32, 2, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
x, y = 2, 2

# Move sequence
moves = 'duldrl'

# Function to perform a move
def move(x, y, direction):
    if direction == 'u':
        if x > 0:
            x -= 1
    elif direction == 'd':
        if x < 4:
            x += 1
    elif direction == 'l':
        if y > 0:
            y -= 1
    elif direction == 'r':
        if y < 4:
            y += 1
    return x, y

# Simulate the moves
for move_direction in moves:
    new_x, new_y = move(x, y, move_direction)
    if grid[new_x][new_y] == grid[x][y]:
        grid[new_x][new_y] *= 2
        grid[x][y] = 0
    elif grid[new_x][new_y] == 0:
        grid[new_x][new_y] = grid[x][y]
        grid[x][y] = 0
    x, y = new_x, new_y

# Print the final grid
print(grid)