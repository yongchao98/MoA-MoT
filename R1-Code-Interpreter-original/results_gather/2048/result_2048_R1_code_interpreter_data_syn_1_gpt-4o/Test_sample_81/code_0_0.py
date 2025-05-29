# Initial grid
grid = [
    [0, 0, 0, 0, 2],
    [0, 0, 0, 0, 0],
    [0, 0, 8, 0, 32],
    [0, 0, 0, 0, 4],
    [4, 32, 0, 0, 0]
]

# Initial position of the number to move
x, y = 2, 2

# Movement sequence
moves = 'rdlluu'

# Function to perform a move
def move(x, y, direction):
    if direction == 'r':
        return x, y + 1
    elif direction == 'd':
        return x + 1, y
    elif direction == 'l':
        return x, y - 1
    elif direction == 'u':
        return x - 1, y
    return x, y

# Simulate the movement
for move_direction in moves:
    new_x, new_y = move(x, y, move_direction)
    if grid[new_x][new_y] == grid[x][y]:
        # Combine numbers
        grid[new_x][new_y] *= 2
        grid[x][y] = 0
    elif grid[new_x][new_y] == 0:
        # Move to the new position
        grid[new_x][new_y] = grid[x][y]
        grid[x][y] = 0
    # Update position
    x, y = new_x, new_y

# Print the final grid
print(grid)