# Initial grid
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 2, 8, 0],
    [0, 0, 0, 8, 0],
    [32, 0, 0, 0, 0],
    [0, 8, 0, 0, 0]
]

# Initial position of the number
x, y = 2, 4

# Sequence of moves
moves = 'rlulur'

# Function to perform the moves
def move_number(grid, x, y, moves):
    for move in moves:
        prev_x, prev_y = x, y  # Store previous position

        if move == 'r':  # Move right
            if y < 4:
                y += 1
        elif move == 'l':  # Move left
            if y > 0:
                y -= 1
        elif move == 'u':  # Move up
            if x > 0:
                x -= 1
        elif move == 'd':  # Move down
            if x < 4:
                x += 1
        
        # Check for combination
        if grid[x][y] == grid[prev_x][prev_y]:
            grid[x][y] *= 2
            grid[prev_x][prev_y] = 0

    return grid

# Update the grid
updated_grid = move_number(grid, x, y, moves)

# Print the updated grid
print(updated_grid)