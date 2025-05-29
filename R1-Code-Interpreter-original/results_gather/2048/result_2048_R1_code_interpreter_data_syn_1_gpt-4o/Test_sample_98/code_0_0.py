# Initial grid
grid = [
    [0, 0, 0, 0, 0],
    [4, 0, 0, 0, 0],
    [0, 0, 0, 0, 2],
    [0, 0, 16, 0, 4],
    [0, 0, 8, 16, 16]
]

# Initial position of the number to move
position = (4, 5)

# Move sequence
moves = 'rudlrd'

# Function to perform a move
def move_number(grid, position, direction):
    x, y = position
    if direction == 'r' and y < 4:
        y += 1
    elif direction == 'l' and y > 0:
        y -= 1
    elif direction == 'u' and x > 0:
        x -= 1
    elif direction == 'd' and x < 4:
        x += 1
    return (x, y)

# Simulate the moves
for move in moves:
    new_position = move_number(grid, position, move)
    x, y = new_position
    px, py = position
    if grid[x][y] == grid[px][py]:
        grid[x][y] *= 2
        grid[px][py] = 0
    elif grid[x][y] == 0:
        grid[x][y] = grid[px][py]
        grid[px][py] = 0
    position = new_position

# Print the final grid
print(grid)