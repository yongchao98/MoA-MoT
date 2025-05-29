# Initial grid
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 2, 8, 0],
    [0, 0, 0, 8, 0],
    [32, 0, 0, 0, 0],
    [0, 8, 0, 0, 0]
]

# Initial position of the number to move
position = (2, 3)  # Adjusted to 0-based index

# Sequence of moves
moves = 'rlulur'

# Function to perform a move
def move_and_combine(grid, position, direction):
    x, y = position
    if direction == 'r':
        while y < 4:
            if grid[x][y+1] == 0:
                grid[x][y+1] = grid[x][y]
                grid[x][y] = 0
                y += 1
            elif grid[x][y+1] == grid[x][y]:
                grid[x][y+1] *= 2
                grid[x][y] = 0
                y += 1
                break
            else:
                break
    elif direction == 'l':
        while y > 0:
            if grid[x][y-1] == 0:
                grid[x][y-1] = grid[x][y]
                grid[x][y] = 0
                y -= 1
            elif grid[x][y-1] == grid[x][y]:
                grid[x][y-1] *= 2
                grid[x][y] = 0
                y -= 1
                break
            else:
                break
    elif direction == 'u':
        while x > 0:
            if grid[x-1][y] == 0:
                grid[x-1][y] = grid[x][y]
                grid[x][y] = 0
                x -= 1
            elif grid[x-1][y] == grid[x][y]:
                grid[x-1][y] *= 2
                grid[x][y] = 0
                x -= 1
                break
            else:
                break
    elif direction == 'd':
        while x < 4:
            if grid[x+1][y] == 0:
                grid[x+1][y] = grid[x][y]
                grid[x][y] = 0
                x += 1
            elif grid[x+1][y] == grid[x][y]:
                grid[x+1][y] *= 2
                grid[x][y] = 0
                x += 1
                break
            else:
                break
    return (x, y)

# Process each move
for move in moves:
    position = move_and_combine(grid, position, move)

# Print the final grid
print(grid)