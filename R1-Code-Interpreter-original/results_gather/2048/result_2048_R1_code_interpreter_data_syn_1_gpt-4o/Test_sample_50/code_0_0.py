# Initial matrix
matrix = [
    [4, 8, 0, 8],
    [16, 0, 0, 0],
    [16, 2, 0, 0],
    [0, 32, 16, 0]
]

# Initial position of the number to move
position = (3, 0)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'udul'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'u':
        return (max(row - 1, 0), col)
    elif direction == 'd':
        return (min(row + 1, 3), col)
    elif direction == 'l':
        return (row, max(col - 1, 0))
    elif direction == 'r':
        return (row, min(col + 1, 3))
    return position

# Simulate the moves
for move_direction in moves:
    new_position = move(position, move_direction)
    # Check if we can combine
    if matrix[new_position[0]][new_position[1]] == matrix[position[0]][position[1]]:
        matrix[new_position[0]][new_position[1]] *= 2
        matrix[position[0]][position[1]] = 0
    else:
        # Move the number
        matrix[new_position[0]][new_position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
    position = new_position

# Print the final matrix
print(matrix)