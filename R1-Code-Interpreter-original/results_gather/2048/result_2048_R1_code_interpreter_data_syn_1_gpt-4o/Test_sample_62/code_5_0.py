# Initial grid
grid = [
    [8, 0, 0, 32, 0],
    [4, 0, 0, 16, 0],
    [0, 0, 0, 0, 0],
    [2, 4, 0, 2, 4],
    [4, 0, 0, 0, 0]
]

# Correct initial position of the number to move (0-indexed)
position = (3, 3)

# Sequence of moves
moves = 'dldur'

# Function to perform a move
def move(grid, position, direction):
    x, y = position
    if direction == 'u':
        new_x, new_y = x - 1, y
    elif direction == 'd':
        new_x, new_y = x + 1, y
    elif direction == 'l':
        new_x, new_y = x, y - 1
    elif direction == 'r':
        new_x, new_y = x, y + 1
    else:
        return position  # Invalid direction, no move

    # Check if the new position is within bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        # Check if we can combine numbers
        if grid[new_x][new_y] == grid[x][y]:
            grid[new_x][new_y] *= 2
            grid[x][y] = 0
        elif grid[new_x][new_y] == 0:
            grid[new_x][new_y] = grid[x][y]
            grid[x][y] = 0
        return (new_x, new_y)
    return position

# Perform the sequence of moves
for move_direction in moves:
    position = move(grid, position, move_direction)

# Print the final grid
print(grid)