# Initial grid setup
grid = [
    [32, 0, 32, 8, 0],
    [16, 0, 0, 0, 0],
    [0, 0, 16, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 2, 0, 0, 0]
]

# Initial position of the number
position = (4, 1)  # (row, column) in 0-indexed

# Move sequence
moves = 'llldud'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'u' and row > 0:
        return (row - 1, col)
    elif direction == 'd' and row < 4:
        return (row + 1, col)
    elif direction == 'l' and col > 0:
        return (row, col - 1)
    elif direction == 'r' and col < 4:
        return (row, col + 1)
    return position

# Simulate the moves
for direction in moves:
    new_position = move(position, direction)
    # Check if we can combine
    if grid[new_position[0]][new_position[1]] == grid[position[0]][position[1]]:
        grid[new_position[0]][new_position[1]] *= 2
        grid[position[0]][position[1]] = 0
    elif grid[new_position[0]][new_position[1]] == 0:
        grid[new_position[0]][new_position[1]] = grid[position[0]][position[1]]
        grid[position[0]][position[1]] = 0
    position = new_position

# Print the final grid
print(grid)