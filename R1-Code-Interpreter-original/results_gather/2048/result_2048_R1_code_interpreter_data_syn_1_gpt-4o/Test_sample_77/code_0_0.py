# Initialize the grid
grid = [
    [0, 0, 0, 8, 0],
    [0, 4, 0, 0, 0],
    [0, 0, 2, 0, 8],
    [0, 0, 0, 0, 0],
    [0, 16, 0, 0, 0]
]

# Starting position of the number 16
position = (4, 1)  # (row, column) in 0-indexed

# Move sequence
moves = 'lrrul'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'u':
        return (max(row - 1, 0), col)
    elif direction == 'd':
        return (min(row + 1, 4), col)
    elif direction == 'l':
        return (row, max(col - 1, 0))
    elif direction == 'r':
        return (row, min(col + 1, 4))
    return position

# Simulate the moves
for move_direction in moves:
    new_position = move(position, move_direction)
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