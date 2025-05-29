# Initial matrix
matrix = [
    [0, 0, 0, 0, 16],
    [4, 32, 0, 32, 4],
    [2, 0, 0, 0, 0],
    [0, 0, 0, 2, 8],
    [0, 0, 0, 8, 0]
]

# Initial position of the number
position = (0, 4)

# Sequence of moves
moves = 'luuddd'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    value = matrix[x][y]
    matrix[x][y] = 0  # Clear the initial position

    for move in moves:
        if move == 'l' and y > 0:
            y -= 1
        elif move == 'r' and y < 4:
            y += 1
        elif move == 'u' and x > 0:
            x -= 1
        elif move == 'd' and x < 4:
            x += 1

        # Check for merging
        if matrix[x][y] == value:
            value *= 2
            matrix[x][y] = 0  # Clear the merged position
        elif matrix[x][y] == 0:
            matrix[x][y] = value
        else:
            # If the position is occupied by a different number, stop moving
            break

    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)