# Initial grid
grid = [
    [0, 4, 8, 0],
    [0, 2, 0, 0],
    [0, 32, 2, 0],
    [2, 4, 0, 0]
]

# Initial position of the number to move
position = (2, 2)

# Sequence of moves
moves = 'rrud'

# Function to perform the moves
def move_number(grid, position, moves):
    x, y = position
    for move in moves:
        # Clear the current position
        grid[x][y] = 0
        if move == 'r':
            if y < 3:  # Move right
                y += 1
        elif move == 'l':
            if y > 0:  # Move left
                y -= 1
        elif move == 'u':
            if x > 0:  # Move up
                x -= 1
        elif move == 'd':
            if x < 3:  # Move down
                x += 1
        # Place the number in the new position
        grid[x][y] = 32
    return grid

# Perform the moves
final_grid = move_number(grid, position, moves)

# Print the final grid
print(final_grid)