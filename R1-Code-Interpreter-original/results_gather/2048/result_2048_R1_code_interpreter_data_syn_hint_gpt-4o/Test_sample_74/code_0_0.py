# Initial matrix
matrix = [
    [0, 0, 0, 0, 0],
    [0, 0, 2, 0, 32],
    [0, 0, 0, 0, 0],
    [0, 0, 16, 0, 0],
    [16, 2, 0, 0, 0]
]

# Initial position of the number to move
position = (1, 4)  # (row, column) in 0-indexed
value = matrix[position[0]][position[1]]

# Sequence of moves
moves = 'rrurl'

# Function to move the number in the matrix
def move_number(matrix, position, value, moves):
    for move in moves:
        if move == 'r':
            # Move right
            if position[1] < 4:
                if matrix[position[0]][position[1] + 1] == value:
                    # Combine numbers
                    matrix[position[0]][position[1] + 1] *= 2
                    value = matrix[position[0]][position[1] + 1]
                    matrix[position[0]][position[1]] = 0
                elif matrix[position[0]][position[1] + 1] == 0:
                    # Move to empty space
                    matrix[position[0]][position[1] + 1] = value
                    matrix[position[0]][position[1]] = 0
                position = (position[0], position[1] + 1)
        elif move == 'l':
            # Move left
            if position[1] > 0:
                if matrix[position[0]][position[1] - 1] == value:
                    # Combine numbers
                    matrix[position[0]][position[1] - 1] *= 2
                    value = matrix[position[0]][position[1] - 1]
                    matrix[position[0]][position[1]] = 0
                elif matrix[position[0]][position[1] - 1] == 0:
                    # Move to empty space
                    matrix[position[0]][position[1] - 1] = value
                    matrix[position[0]][position[1]] = 0
                position = (position[0], position[1] - 1)
        elif move == 'u':
            # Move up
            if position[0] > 0:
                if matrix[position[0] - 1][position[1]] == value:
                    # Combine numbers
                    matrix[position[0] - 1][position[1]] *= 2
                    value = matrix[position[0] - 1][position[1]]
                    matrix[position[0]][position[1]] = 0
                elif matrix[position[0] - 1][position[1]] == 0:
                    # Move to empty space
                    matrix[position[0] - 1][position[1]] = value
                    matrix[position[0]][position[1]] = 0
                position = (position[0] - 1, position[1])
        elif move == 'd':
            # Move down
            if position[0] < 4:
                if matrix[position[0] + 1][position[1]] == value:
                    # Combine numbers
                    matrix[position[0] + 1][position[1]] *= 2
                    value = matrix[position[0] + 1][position[1]]
                    matrix[position[0]][position[1]] = 0
                elif matrix[position[0] + 1][position[1]] == 0:
                    # Move to empty space
                    matrix[position[0] + 1][position[1]] = value
                    matrix[position[0]][position[1]] = 0
                position = (position[0] + 1, position[1])
    return matrix

# Move the number according to the sequence
result_matrix = move_number(matrix, position, value, moves)

# Print the resulting matrix
print(result_matrix)