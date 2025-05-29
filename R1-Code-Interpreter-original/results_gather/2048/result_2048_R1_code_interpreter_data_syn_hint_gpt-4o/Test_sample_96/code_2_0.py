# Initial matrix
matrix = [
    [0, 0, 0, 0, 32],
    [0, 16, 0, 0, 16],
    [2, 0, 8, 0, 0],
    [4, 0, 0, 0, 4],
    [0, 0, 8, 0, 4]
]

# Initial position of the number to move
position = (4, 2)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'lululr'

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
    return position

# Perform the moves
for move_direction in moves:
    new_position = move(position, move_direction)
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        # Check if the numbers can be combined
        if matrix[new_position[0]][new_position[1]] == matrix[position[0]][position[1]]:
            # Combine the numbers
            matrix[new_position[0]][new_position[1]] *= 2
            matrix[position[0]][position[1]] = 0
            # Update the position
            position = new_position
        elif matrix[new_position[0]][new_position[1]] == 0:
            # Move the number
            matrix[new_position[0]][new_position[1]] = matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
            # Update the position
            position = new_position

# Print the final matrix
for row in matrix:
    print(row)