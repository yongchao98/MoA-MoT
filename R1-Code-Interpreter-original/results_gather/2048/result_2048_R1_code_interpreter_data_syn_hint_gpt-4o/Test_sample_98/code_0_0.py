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

# Sequence of moves
moves = 'rudlrd'

# Function to perform a move
def move_number(grid, position, direction):
    x, y = position
    if direction == 'r':
        new_x, new_y = x, y + 1
    elif direction == 'u':
        new_x, new_y = x - 1, y
    elif direction == 'd':
        new_x, new_y = x + 1, y
    elif direction == 'l':
        new_x, new_y = x, y - 1
    else:
        return position  # Invalid direction, return the same position

    # Check if the new position is within bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        # Check if the numbers can combine
        if grid[new_x][new_y] == grid[x][y]:
            grid[new_x][new_y] *= 2
            grid[x][y] = 0
        elif grid[new_x][new_y] == 0:
            grid[new_x][new_y] = grid[x][y]
            grid[x][y] = 0
        return (new_x, new_y)
    return position

# Perform the sequence of moves
for move in moves:
    position = move_number(grid, position, move)

# Print the final grid
print(grid)