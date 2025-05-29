# Initial matrix
matrix = [
    [0, 0, 0, 32, 0],
    [0, 0, 0, 0, 16],
    [0, 2, 16, 0, 32],
    [0, 0, 2, 2, 0],
    [0, 32, 16, 0, 0]
]

# Initial position of the number to move
position = (2, 4)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'rulrdr'

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
for move_direction in moves:
    new_position = move(position, move_direction)
    # Check if we can combine
    if matrix[new_position[0]][new_position[1]] == matrix[position[0]][position[1]]:
        matrix[new_position[0]][new_position[1]] *= 2
        matrix[position[0]][position[1]] = 0
    elif matrix[new_position[0]][new_position[1]] == 0:
        matrix[new_position[0]][new_position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
    position = new_position

# Print the final matrix
print(matrix)