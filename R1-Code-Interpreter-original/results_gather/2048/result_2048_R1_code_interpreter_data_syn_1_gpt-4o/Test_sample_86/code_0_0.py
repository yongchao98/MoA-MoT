# Initial grid
grid = [
    [0, 16, 0, 0, 32],
    [0, 0, 0, 4, 0],
    [0, 0, 0, 0, 32],
    [8, 0, 0, 0, 0],
    [0, 32, 4, 8, 0]
]

# Initial position of the number to move
position = (2, 4)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'dludrd'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'u':
        return (row - 1, col)
    elif direction == 'd':
        return (row + 1, col)
    elif direction == 'l':
        return (row, col - 1)
    elif direction == 'r':
        return (row, col + 1)

# Perform the moves
for direction in moves:
    new_position = move(position, direction)
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        # Check if the numbers can combine
        if grid[new_position[0]][new_position[1]] == grid[position[0]][position[1]]:
            grid[new_position[0]][new_position[1]] *= 2
            grid[position[0]][position[1]] = 0
        elif grid[new_position[0]][new_position[1]] == 0:
            grid[new_position[0]][new_position[1]] = grid[position[0]][position[1]]
            grid[position[0]][position[1]] = 0
        position = new_position

# Print the final grid
print(grid)