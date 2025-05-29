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
    return x, y

# Perform the moves
final_position = move_number(grid, position, moves)

# Update the grid
grid[2][2] = 0  # Clear the original position
grid[final_position[0]][final_position[1]] = 32  # Place the number in the new position

# Print the final grid
print(grid)