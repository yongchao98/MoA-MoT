# Initial matrix
matrix = [
    [0, 4, 0, 4],
    [0, 0, 0, 0],
    [0, 4, 16, 0],
    [0, 0, 0, 0]
]

# Function to move the matrix in a given direction
def move(matrix, direction):
    if direction == 'd':
        # Move down
        for col in range(4):
            for row in range(2, -1, -1):
                if matrix[row][col] != 0:
                    for next_row in range(row + 1, 4):
                        if matrix[next_row][col] == 0:
                            matrix[next_row][col] = matrix[next_row - 1][col]
                            matrix[next_row - 1][col] = 0
                        elif matrix[next_row][col] == matrix[next_row - 1][col]:
                            matrix[next_row][col] *= 2
                            matrix[next_row - 1][col] = 0
                            break
                        else:
                            break
    elif direction == 'u':
        # Move up
        for col in range(4):
            for row in range(1, 4):
                if matrix[row][col] != 0:
                    for next_row in range(row - 1, -1, -1):
                        if matrix[next_row][col] == 0:
                            matrix[next_row][col] = matrix[next_row + 1][col]
                            matrix[next_row + 1][col] = 0
                        elif matrix[next_row][col] == matrix[next_row + 1][col]:
                            matrix[next_row][col] *= 2
                            matrix[next_row + 1][col] = 0
                            break
                        else:
                            break

# Sequence of moves
sequence = 'duuu'

# Apply each move in the sequence
for move_direction in sequence:
    move(matrix, move_direction)

# Print the final matrix
print(matrix)