# Initial matrix
matrix = [
    [0, 16, 0, 0, 32],
    [0, 0, 0, 4, 0],
    [0, 0, 0, 0, 32],
    [8, 0, 0, 0, 0],
    [0, 32, 4, 8, 0]
]

# Initial position of the number to move
position = (2, 4)  # 0-indexed

# Sequence of moves
moves = 'dludrd'

# Function to perform a move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'u':
        new_x, new_y = x - 1, y
    elif direction == 'd':
        new_x, new_y = x + 1, y
    elif direction == 'l':
        new_x, new_y = x, y - 1
    elif direction == 'r':
        new_x, new_y = x, y + 1
    else:
        return position  # Invalid direction, no move

    # Check if the new position is within bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        # Check if the numbers can combine
        if matrix[new_x][new_y] == matrix[x][y]:
            matrix[new_x][new_y] *= 2
            matrix[x][y] = 0
        elif matrix[new_x][new_y] == 0:
            matrix[new_x][new_y] = matrix[x][y]
            matrix[x][y] = 0
        return (new_x, new_y)
    return position

# Perform the sequence of moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)