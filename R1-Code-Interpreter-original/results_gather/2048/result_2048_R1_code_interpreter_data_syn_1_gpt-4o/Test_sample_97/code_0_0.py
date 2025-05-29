# Initial grid
grid = [
    [0, 4, 0, 2, 0],
    [0, 16, 8, 0, 32],
    [0, 0, 0, 16, 0],
    [4, 16, 0, 0, 0],
    [0, 8, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 4)

# Movement sequence
moves = 'ulurrl'

# Function to perform a move
def move_number(grid, position, direction):
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
        # Check if the numbers can be combined
        if grid[new_x][new_y] == grid[x][y]:
            grid[new_x][new_y] *= 2
            grid[x][y] = 0
        else:
            # Move the number to the new position
            grid[new_x][new_y], grid[x][y] = grid[x][y], grid[new_x][new_y]
        return new_x, new_y
    return position

# Simulate the moves
for move in moves:
    position = move_number(grid, position, move)

# Print the final grid
print(grid)